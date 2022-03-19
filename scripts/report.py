import json
import shutil
import pandas as pd
import mols2grid
from utils import ROOT, DATA, RESULTS, N_WORKERS


RESULTS.mkdir(parents=True, exist_ok=True)

with open(DATA / ".fetched_count") as fi_fetched, \
     open(DATA / ".processed_unique_count") as fi_processed, \
     open(DATA / ".failed_count") as fi_failed, \
     open(DATA / ".timestamp") as fi_timestamp, \
     open(DATA / ".timing") as fi_timing, \
     open(RESULTS / "results.json", "w") as fo:

    n_mols = int(fi_processed.read())
    n_fails = int(fi_failed.read())
    acc = 100 * (n_mols - n_fails) / n_mols
    results = {
        "Date": fi_timestamp.read(),
        "Accuracy (%)": acc,
        "Number of molecules fetched": int(fi_fetched.read()),
        "Number of molecules processed": n_mols,
        "Number of molecules failed": n_fails,
        "Timing (s)": float(fi_timing.read()),
        "Number of threads": N_WORKERS,
    }
    json.dump(results, fo)

# copy failures to results folder
shutil.copy(DATA / "chembl_failed.smi", RESULTS / "chembl_failed.smi")

# generate grid of failed molecules
in_file = DATA / "chembl_failed.smi"
df = pd.read_csv(in_file, sep=" ", names=["SMILES", "ChEMBL id"])

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
    df, output=RESULTS / "failed_molecules.html",
    subset=["ChEMBL id", "img"],
    tooltip=["SMILES"], tooltip_trigger="hover",
    clearBackground=False,
    callback=callback,
)

print(f"Wrote reports in {RESULTS.relative_to(ROOT)}")
