import json
import shutil

import mols2grid
import pandas as pd
from MDAnalysis import __version__ as mda_version
from utils import DATA, RESULTS, ROOT


def make_report():
    with open(DATA / ".fetched_count") as fi_fetched, open(
        DATA / ".processed_unique_count"
    ) as fi_processed, open(DATA / ".failed_count") as fi_failed, open(
        DATA / ".timestamp"
    ) as fi_timestamp, open(
        DATA / ".timing"
    ) as fi_timing, open(
        RESULTS / "results.json", "w"
    ) as fo:
        n_mols = int(fi_processed.read())
        n_fails = int(fi_failed.read())
        accuracy = 100 * (n_mols - n_fails) / n_mols
        results = {
            "Date": fi_timestamp.read(),
            "MDAnalysis version": mda_version,
            "Accuracy (%)": accuracy,
            "Number of molecules fetched": int(fi_fetched.read()),
            "Number of molecules processed": n_mols,
            "Number of molecules failed": n_fails,
            "Timing (s)": float(fi_timing.read()),
        }
        json.dump(results, fo, indent=4)

    # copy failures to results folder
    shutil.copy(DATA / "chembl_failed.smi", RESULTS / "failed_molecules.smi")
    return accuracy


if __name__ == "__main__":
    print("Making reports")
    RESULTS.mkdir(parents=True, exist_ok=True)

    accuracy = make_report()

    # generate grid of failed molecules
    in_file = DATA / "chembl_failed.smi"
    df = pd.read_csv(in_file, sep=" ", names=["SMILES", "ChEMBL id"])

    callback = mols2grid.make_popup_callback(
        title="${data['ChEMBL-id']}",
        js="""
            var mol = RDKitModule.get_mol(data["SMILES"]);
            var svg = mol.get_svg(600, 500);
            mol.delete();
        """,
        html="""
            <div>${svg}</div>
            <p>${data['SMILES']}</p>
        """,
        style="max-width: 620px;",
    )

    mols2grid.save(
        df,
        size=(200, 160),
        n_rows=5,
        n_cols=8,
        output=RESULTS / "failed_molecules.html",
        subset=["ChEMBL id", "img"],
        tooltip=["SMILES"],
        tooltip_trigger="hover",
        clearBackground=False,
        callback=callback,
    )

    # generate badge
    with open(RESULTS / "badge.json", "w") as f:
        badge = {
            "schemaVersion": 1,
            "label": "accuracy",
            "message": f"{accuracy:.2f}%",
            "color": "success",
        }
        json.dump(badge, f)

    print(f"Wrote reports in {RESULTS.relative_to(ROOT)}")
