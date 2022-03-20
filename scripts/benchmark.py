import time
import gzip
import warnings
import multiprocessing as mp
from tqdm.auto import tqdm
from rdkit import Chem, RDLogger
from utils import (ROOT, DATA, N_WORKERS, timestamp,
                   add_Hs_remove_bo_and_charges, enumerate_reordered_mol,
                   assign_bond_orders_and_charges, same_molecules)


# ignore warnings and errors
warnings.filterwarnings("ignore")
RDLogger.DisableLog('rdApp.*')

# files
in_file = DATA / "chembl_processed_unique.smi.gz"
out_file = DATA / "chembl_failed.smi"
num_entries = int(open(DATA / ".processed_unique_count").read()) 
(DATA / ".timestamp").write_text(timestamp)


def validate_entry(smi):
    ref = Chem.MolFromSmiles(smi)
    stripped_mol = add_Hs_remove_bo_and_charges(ref)
    for mol in enumerate_reordered_mol(stripped_mol):
        mol = assign_bond_orders_and_charges(mol)
        if not mol:
            return smi
        mol = Chem.RemoveHs(mol)
        valid = same_molecules(mol, ref)
        if not valid:
            return smi

print("Starting benchmark")
count = 0
start = time.perf_counter()
with mp.Pool(N_WORKERS) as pool, \
     gzip.open(in_file, "rt") as fi, \
     open(out_file, "w") as fo, \
     tqdm(total=num_entries) as pbar:

    for i, result in enumerate(pool.imap_unordered(validate_entry, fi)):
        if result:
            fo.write(result)
            count += 1
            acc = 100 * (num_entries - count) / num_entries
            pbar.set_description(f"Accuracy: {acc:.2f}% | Progress")
        pbar.update()

stop = time.perf_counter()


print(f"Wrote '{out_file.relative_to(ROOT)}' with {count:,} failed entries")
acc = 100 * (num_entries - count) / num_entries
print(f"Final accuracy: {acc:.3f}%")
(DATA / ".failed_count").write_text(str(count))
(DATA / ".timing").write_text(str(stop - start))
