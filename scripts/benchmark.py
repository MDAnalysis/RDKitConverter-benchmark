import gzip
from joblib import Parallel, delayed
from tqdm.auto import tqdm
from rdkit import Chem
from utils import (datapath, timestamp, delimiter, add_Hs_remove_bo_and_charges, 
                   enumerate_reordered_mol, assign_bond_orders_and_charges,
                   is_isomorphic_or_resonance_structure)


mol_file = datapath / f"chembl_{timestamp}_prepared.rdkit.gz"
out_file = datapath / f"chembl_{timestamp}_failed.smi"
delimiter_len = len(delimiter)
num_entries = int(open(datapath / f".molecules_count").read())

def mol_supplier(handle):
    buffer = b""
    for line in handle:
        buffer += line
        if line.endswith(delimiter):
            yield Chem.Mol(buffer[:-delimiter_len])
            buffer = b""    

def validate_mol(mol):
    stripped_mol = add_Hs_remove_bo_and_charges(mol)
    for m in enumerate_reordered_mol(stripped_mol):
        inferred = assign_bond_orders_and_charges(m)
        inferred = Chem.RemoveHs(inferred)
        iso = is_isomorphic_or_resonance_structure(inferred, mol)
        if not iso:
            return False, Chem.MolToSmiles(mol), mol.GetProp("_Name")
    return True

with gzip.open(mol_file, "rb") as fi, \
     open(out_file, "w") as fo:

    count = 0
    for mol in tqdm(mol_supplier(fi), total=num_entries):
        passed = validate_mol(mol)
        if not passed is True:
            _, smi, name = passed
            fo.write(f"{smi} {name}\n")
            count += 1

print(f"Wrote '{out_file}' with {count:,} failed entries")
acc = 100 * (num_entries - count) / num_entries
print(f"Accuracy: {acc:.2f}%")
(datapath / ".failures_count").write_text(str(count))
