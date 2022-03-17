import pandas as pd
import mols2grid
from utils import datapath, timestamp

out_file = datapath / f"chembl_{timestamp}_failed.smi"

# generate grid of failed molecules
df = pd.read_csv(out_file, sep=" ", names=["SMILES", "ChEMBL id"])

callback = mols2grid.make_popup_callback(
    title="${data['ChEMBL-id']}",
    js="""
        var mol = RDKitModule.get_mol(data["SMILES"]);
        var svg = mol.get_svg(600, 600);
        mol.delete();
    """,
    html="""
        <div>${svg}</div>
        <p>${data['SMILES']}</p>
    """,
    style="max-width: 80%;",
)

mols2grid.save(
    df, output=datapath / "failures.html",
    subset=["ChEMBL id", "img"],
    tooltip=["SMILES"], tooltip_trigger="hover",
    clearBackground=False,
    callback=callback,
)
