import gzip
from tqdm.auto import tqdm
from chembl_webresource_client.settings import Settings
Settings.Instance().MAX_LIMIT = 1000
from chembl_webresource_client.new_client import new_client
from utils import datapath, timestamp


print("Querying ChEMBL for molecules with less than 50 heavy atoms")

# fetch molecules from chembl
query = (
    new_client
    .molecule
    .filter(
        molecule_properties__heavy_atoms__lt=50,
    )
    .only([
        "molecule_chembl_id",
        "molecule_structures",
    ])
)
n_entries = len(query)
print(f"Found {n_entries:,} entries")

datapath.mkdir(parents=True, exist_ok=True)
out_file = datapath / f"chembl_{timestamp}_filtered.smi.gz"

with gzip.open(out_file, "wt") as f:
    count = 0
    for entry in tqdm(query, total=n_entries):
        name = entry["molecule_chembl_id"]
        smiles = entry["molecule_structures"]["canonical_smiles"]
        if name and smiles:
            f.write(f"{smiles} {name}\n")
            count += 1

print(f"Wrote '{out_file}' with {count:,} entries")
(datapath / ".chembl_filtered_count").write_text(str(count))
