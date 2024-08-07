import lzma
import multiprocessing as mp
import pickle
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Optional, Sequence, Type, TypeVar

import mols2grid
import networkx as nx
import numpy as np
import pandas as pd
import scaffoldgraph as sg
from benchmark import validate_inferring
from bokeh.io import output_file, save
from bokeh.models import (
    Circle,
    CustomJSHover,
    HoverTool,
    Legend,
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

bokeh_template = """
{% block postamble %}
    <script src="https://unpkg.com/@rdkit/rdkit@2023.3.2-1.0.0/Code/MinimalLib/dist/RDKit_minimal.js"></script>
    <script>
        window
            .initRDKitModule()
            .then(function (RDKit) {
                window.RDKit = RDKit;
                console.log("RDKit version: " + RDKit.version());
            });
    </script>
{% endblock %}"""

js_mol_from_smiles = """
const smiles = value;
const mol = RDKit.get_mol(smiles);
let svg = "";
if (mol.is_valid()) {
    svg = mol.get_svg(%(width)d, %(height)d);
}
mol.delete();
if (svg == "") {
    svg = '<svg width="%(width)d" height="%(height)d" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 %(width)d %(height)d"></svg>';
}
return svg;
"""


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
    def _additional_data(smiles) -> tuple[str, dict]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Something wrong with: {smiles}")
            return smiles, {
                "smiles": smiles,
                "heavy_atom_count": 0,
                "status": "fail",
            }
        return smiles, {
            "smiles": smiles,
            "heavy_atom_count": mol.GetNumHeavyAtoms(),
            "status": "pass" if validate_inferring(mol) else "fail",
        }

    def add_data(self) -> None:
        graph = self.scaffold_network
        added = {}
        with mp.Pool(N_WORKERS) as pool:
            for idx, data in tqdm(
                pool.imap_unordered(self._additional_data, graph.get_scaffold_nodes()),
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
        root = set(failing)
        for failed_idx in failing:
            for _, data in graph.get_parent_scaffolds(failed_idx, data=True):
                if data["status"] == "fail":
                    root.discard(failed_idx)
                    break
        values = {idx: "root" for idx in root}
        nx.set_node_attributes(graph, values, "status")
        with open(smiles_output, "w") as fh:
            fmt_str = "{smiles} {idx} {count}\n"
            for idx in root:
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
        size: tuple[int, int] = (800, 400),
        palette: dict[str, str] = {
            "pass": "#4393C3",
            "fail": "#D6604D",
            "root": "#FA9C4A",
        },
        size_by: str = "mol_count",
        size_bins: Any = "auto",
        size_scale: float = 2,
        node_alpha: float = 0.6,
        edge_alpha: float = 0,
        edge_width: float = 0,
    ) -> None:
        """Interactive bokeh plot of the network"""

        plot = figure(
            tools="pan,tap,wheel_zoom,save,reset",
            active_scroll="wheel_zoom",
            x_axis_location=None,
            y_axis_location=None,
            title=title,
            width=size[0],
            height=size[1],
            sizing_mode="scale_both",
        )
        max_width = f"{self.img_size[0] + 5}px"
        hover_tooltip = (
            """<div style="font-size: 14px; max-width: %s">
        <div style="position: relative; left=0">@smiles{custom}</div>
        <div style="position: relative;">
        <span style="position: relative; font-weight: bold;">@mol_count</span>
        <span style="font-size: 10px; word-wrap: break-word;">@smiles</span>
        </div>
        </div>"""
            % max_width
        )

        plot.add_tools(
            HoverTool(
                tooltips=hover_tooltip,
                formatters={
                    "@smiles": CustomJSHover(
                        code=js_mol_from_smiles
                        % {"width": self.img_size[0], "height": self.img_size[1]}
                    )
                },
            )
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
        min_node_size = 4
        assignments *= size_scale
        assignments += min_node_size
        assignments[indices[-1]] += min_node_size  # benzene slightly larger
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
            size="node_size", fill_color=cmap, line_width=0, line_color="#000000"
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

        legend = Legend(
            title="Status",
            items=[
                (key, [plot.circle(fill_color=color, line_width=0)])
                for key, color in palette.items()
            ],
        )
        plot.add_layout(legend)

        output_file(filename=filename, title=title)
        save(plot, filename=filename, template=bokeh_template)


if __name__ == "__main__":
    suppl = Chem.SmilesMolSupplier(str(RESULTS / "failed_molecules.smi"))
    scaffold_options = ScaffoldOptions()
    # HierS kind to avoid splitting fused rings
    net = Network(kind=sg.HierS, scaffold_options=scaffold_options)
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
        clearBackground=False,
    )
    print("Generating layout and interactive plot")
    net.plot(
        filename=str(RESULTS / "scaffold_network.html"),
        layout=nx.nx_agraph.graphviz_layout,
        layout_kwargs={"prog": "sfdp"},
        size_bins="doane",
    )
