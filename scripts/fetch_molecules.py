import os
import gzip
from tqdm.auto import tqdm
from chembl_webresource_client.settings import Settings
Settings.Instance().MAX_LIMIT = 1000
from chembl_webresource_client.new_client import new_client
from utils import ROOT, DATA, timestamp


# fetch molecules from chembl
min_heavy_atoms = os.getenv("MIN_ATOMS", 2)
max_heavy_atoms = os.getenv("MAX_ATOMS", 50)
print(f"Querying ChEMBL for molecules with between {min_heavy_atoms} and "
      f"{max_heavy_atoms} heavy atoms")

query = (
    new_client
    .molecule
    .filter(
        molecule_properties__heavy_atoms__gte=min_heavy_atoms,
        molecule_properties__heavy_atoms__lte=max_heavy_atoms,
    )
    .only([
        "molecule_chembl_id",
        "molecule_structures",
    ])
)
n_entries = len(query)
print(f"Found {n_entries:,} entries")

DATA.mkdir(parents=True, exist_ok=True)
(DATA / ".timestamp").write_text(timestamp)
out_file = DATA / "chembl_fetched.smi.gz"

with gzip.open(out_file, "wt") as f:
    count = 0
    for entry in tqdm(query, total=n_entries, desc="Fetching ChEMBL"):
        name = entry["molecule_chembl_id"]
        try:
            smiles = entry["molecule_structures"]["canonical_smiles"]
        except KeyError:
            continue
        else:
            f.write(f"{smiles} {name}\n")
            count += 1

print(f"Wrote '{out_file.relative_to(ROOT)}' with {count:,} entries")
(DATA / ".fetched_count").write_text(str(count))
