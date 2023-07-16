import lzma
import multiprocessing as mp
import pickle
from base64 import b64encode
from dataclasses import asdict, dataclass
from functools import partial
from pathlib import Path
from typing import Any, Callable, Optional, Sequence, Type, TypeVar

import mols2grid
import networkx as nx
import numpy as np
import pandas as pd
import scaffoldgraph as sg
from benchmark import validate_inferring
from bokeh.io import output_file, output_notebook, save, show
from bokeh.models import (
    Circle,
    MultiLine,
    NodesAndAdjacentNodes,
    NodesAndLinkedEdges,
    Range1d,
)
from bokeh.plotting import figure, from_networkx
from bokeh.transform import factor_cmap
from rdkit import Chem
from rdkit.Chem import Draw
from tqdm.auto import tqdm
from utils import N_WORKERS, RESULTS


@dataclass()
class ScaffoldOptions:
    """

    Parameters
    ----------
    ring_cutoff : int
        Ignore molecules with more rings than this cutoff to avoid extended
        calculation time.
    progress : bool
        If True display a progress bar to monitor construction progress.
    annotate : bool
        If True write an annotated murcko scaffold SMILES string to each
        molecule edge (molecule --> scaffold).
    flatten_isotopes : bool
        If True remove specific isotopes when initializing the scaffold.
    keep_largest_fragment : bool
        If True when encountering molecules containing disconnected fragments
        initialize the scaffold from only the largest disconnected fragment.
    discharge_and_deradicalize : bool
        If True remove charges and radicals when initializing the scaffold.

    """

    ring_cutoff: int = 10
    progress: bool = True
    annotate: bool = True
    flatten_isotopes: bool = True
    keep_largest_fragment: bool = False
    discharge_and_deradicalize: bool = False

    @property
    def parameters(self) -> dict:
        return asdict(self)


ScaffoldNetworkT = TypeVar("ScaffoldNetworkT", bound=sg.core.ScaffoldGraph)


@dataclass()
class Network:
    kind: Type[ScaffoldNetworkT]
    scaffold_options: ScaffoldOptions
    scaffold_network: ScaffoldNetworkT = None
    use_svg: bool = True
    img_size: tuple[int, int] = (160, 120)
    draw_options: Draw.MolDrawOptions = Draw.MolDrawOptions()

    def _calculate_network(self, mols: Sequence[Chem.Mol]) -> None:
        graph = sg.HierS.from_supplier(mols, **self.scaffold_options.parameters)
        print(
            f"Generated scaffold network from {graph.num_molecule_nodes} molecules "
            f"with {graph.num_scaffold_nodes} scaffolds"
        )
        self.scaffold_network = graph

    @staticmethod
    def draw_svg(mol, img_size: tuple[int, int]) -> str:
        d2d = Draw.MolDraw2DSVG(*img_size)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        return d2d.GetDrawingText().replace("\n", "")

    @staticmethod
    def draw_png(mol, img_size: tuple[int, int]) -> str:
        d2d = Draw.MolDraw2DCairo(*img_size)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        png = d2d.GetDrawingText()
        return b64encode(png).decode("utf-8")

    @staticmethod
    def _additional_data(args) -> tuple[str, dict]:
        smiles, get_img = args
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Something wrong with: {smiles}")
            return smiles, {
                "smiles": smiles,
                "heavy_atom_count": 0,
                "status": "fail",
                "mol_img": get_img(Chem.MolFromSmiles("")),
            }
        return smiles, {
            "smiles": smiles,
            "heavy_atom_count": mol.GetNumHeavyAtoms(),
            "status": "pass" if validate_inferring(mol) else "fail",
            "mol_img": get_img(mol),
        }

    def add_data(self) -> None:
        graph = self.scaffold_network
        added = {}
        img_func = self.draw_svg if self.use_svg else self.draw_png
        get_img = partial(img_func, img_size=self.img_size)
        arguments = ((smiles, get_img) for smiles in graph.get_scaffold_nodes())
        with mp.Pool(N_WORKERS) as pool:
            for idx, data in tqdm(
                pool.imap_unordered(self._additional_data, arguments),
                total=graph.num_scaffold_nodes,
                desc="Calculating additional data",
            ):
                added[idx] = data
        nx.set_node_attributes(graph, added)
        added = {
            smiles: {
                "mol_count": (count := len(graph.get_molecules_for_scaffold(smiles))),
                "weight": count,
            }
            for smiles in graph.get_scaffold_nodes()
        }
        nx.set_node_attributes(graph, added)
        graph.remove_nodes_from([idx for idx in graph.get_molecule_nodes()])

    def generate_network(
        self, mols: Sequence[Chem.Mol], pkl_file: Path, overwrite: bool = False
    ) -> None:
        """Generate and save or load scaffold network"""
        if not overwrite and pkl_file.exists():
            with lzma.open(pkl_file, "rb") as fh:
                graph: ScaffoldNetworkT = pickle.load(fh)
                print(
                    f"Loaded scaffold network with {graph.num_scaffold_nodes} scaffolds"
                )
                self.scaffold_network = graph
                return graph

        self._calculate_network(mols)
        self.add_data()
        nx.relabel_nodes(
            self.scaffold_network,
            {name: idx for idx, name in enumerate(self.scaffold_network.nodes)},
            copy=False,
        )

        with lzma.open(pkl_file, "wb") as fh:
            pickle.dump(self.scaffold_network, fh)

    def extract_failures(self, smiles_output: str) -> None:
        """Extract failing scaffolds that don't have any parent scaffold failing, i.e.
        likely cause of failures"""
        graph = net.scaffold_network
        failing = [
            idx for idx, data in graph.nodes(data=True) if data["status"] == "fail"
        ]
        highlight = set(failing)
        for failed_idx in failing:
            for _, data in graph.get_parent_scaffolds(failed_idx, data=True):
                if data["status"] == "fail":
                    highlight.discard(failed_idx)
        values = {idx: "highlight" for idx in highlight}
        nx.set_node_attributes(graph, values, "status")
        with open(smiles_output, "w") as fh:
            fmt_str = "{smiles} {idx} {count}\n"
            for idx in highlight:
                data = graph.nodes[idx]
                fh.write(
                    fmt_str.format(
                        smiles=data["smiles"], idx=idx, count=data["mol_count"]
                    )
                )

    def plot(
        self,
        *,
        layout: Callable,
        layout_kwargs: dict = {},
        title: str = "Scaffold Network",
        filename: Optional[str] = None,
        size: tuple[int, int] = (1000, 600),
        palette: dict[str, str] = {
            "pass": "#4393C3",
            "fail": "#D6604D",
            "highlight": "#FA9C4A",
        },
        size_by: str = "mol_count",
        size_bins: Any = "auto",
        size_scale: float = 2,
        node_alpha: float = 0.8,
        edge_alpha: float = 0.3,
        edge_width: float = 1,
    ) -> None:
        """Interactive bokeh plot of the network"""
        img_tooltip = (
            "@mol_img{safe}"
            if self.use_svg
            else '<img src="data:image/png;base64,@mol_img{safe}">'
        )
        hover_tooltip = f"""<div style="font-size: 14px; position: relative;">
        <div style="left=0">{img_tooltip}</div>
        <div><span style="font-weight: bold;">@mol_count</span>
        <span style="font-size: 10px;">@smiles</span>
        </div></div>"""

        plot = figure(
            tooltips=hover_tooltip,
            tools="pan,tap,wheel_zoom,save,reset",
            active_scroll="wheel_zoom",
            x_axis_location=None,
            y_axis_location=None,
            title=title,
            width=size[0],
            height=size[1],
        )
        plot.grid.grid_line_color = None
        graph = self.scaffold_network

        # bin data to calculate size of nodes
        values = np.array([data[size_by] for _, data in graph.nodes(data=True)])
        # limit benzene size to second biggest node
        indices = np.argsort(values)
        values[indices[-1]] = values[indices[-2]]
        # make bins and assign
        _, bins = np.histogram(values, bins=size_bins)
        assignments = np.digitize(values, bins, right=True).astype(float)
        # adjust sizes
        min_node_size = 3
        assignments[indices[-1]] += min_node_size  # benzene slightly larger
        assignments *= size_scale
        assignments[assignments == 0] = min_node_size
        size = {idx: s for idx, s in zip(graph.nodes, assignments)}
        nx.set_node_attributes(graph, size, "node_size")

        graph_model = from_networkx(graph, layout, **layout_kwargs)
        coords = np.array(
            list(graph_model.layout_provider.graph_layout.values()), dtype=np.float32
        )
        x_min, y_min = coords.min(axis=0)
        x_max, y_max = coords.max(axis=0)
        plot.x_range = Range1d(x_min, x_max)
        plot.y_range = Range1d(y_min, y_max)
        plot.renderers.append(graph_model)

        cmap = factor_cmap(
            "status", palette=tuple(palette.values()), factors=tuple(palette)
        )

        graph_model.node_renderer.glyph = Circle(
            size="node_size", fill_color=cmap, line_width=0, fill_alpha=node_alpha
        )
        node_interaction_glyph = Circle(
            size="node_size", fill_color=cmap, line_width=2, line_color="#000000"
        )
        graph_model.node_renderer.selection_glyph = node_interaction_glyph
        graph_model.node_renderer.hover_glyph = node_interaction_glyph

        graph_model.edge_renderer.glyph = MultiLine(
            line_alpha=edge_alpha, line_width=edge_width
        )
        edge_interaction_glyph = MultiLine(
            line_alpha=1, line_width=2, line_color="#000000"
        )
        graph_model.edge_renderer.selection_glyph = edge_interaction_glyph
        graph_model.edge_renderer.hover_glyph = edge_interaction_glyph

        graph_model.selection_policy = NodesAndAdjacentNodes()
        graph_model.inspection_policy = NodesAndLinkedEdges()

        if filename:
            output_file(filename=filename, title=title)
            save(plot, filename=filename)
        else:
            output_notebook()
            show(plot)


if __name__ == "__main__":
    suppl = Chem.SmilesMolSupplier(str(RESULTS / "failed_molecules.smi"))
    scaffold_options = ScaffoldOptions()
    # HierS kind to avoid splitting fused rings
    net = Network(kind=sg.HierS, scaffold_options=scaffold_options, use_svg=False)
    net.generate_network(suppl, RESULTS / "scaffold_networkx.pkl.xz", overwrite=False)
    net.extract_failures(str(RESULTS / "failed_scaffolds.smi"))
    df = pd.read_csv(
        str(RESULTS / "failed_scaffolds.smi"),
        names=["SMILES", "idx", "mol count"],
        sep=" ",
    )
    mols2grid.save(
        df,
        subset=["idx", "img"],
        tooltip=["mol count", "SMILES"],
        size=(200, 160),
        n_rows=5,
        n_cols=8,
        sort_by="mol count",
        transform={"mol count": lambda x: -x},
        output=str(RESULTS / "failed_scaffolds.html"),
        tooltip_trigger="hover",
        clearBackground=False,
    )
    print("Generating layout and interactive plot")
    net.plot(
        filename=str(RESULTS / "scaffold_network.html"),
        layout=nx.nx_agraph.graphviz_layout,
        layout_kwargs={"prog": "sfdp"},
        size_bins="doane",
    )
