import gzip
import multiprocessing as mp
import os

from rdkit import Chem, RDLogger
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from tqdm.auto import tqdm

from utils import DATA, N_WORKERS, ROOT, apply_reaction

# ignore warnings and errors
RDLogger.DisableLog("rdApp.*")

# filters
min_heavy_atoms = int(os.getenv("MIN_ATOMS", 2))
max_heavy_atoms = int(os.getenv("MAX_ATOMS", 50))

# sanitization reactions
REACTIONS = [
    "[#6:1]-[N+:2]#[N:3]=[N-:4]>>[#6:1]-[N+0:2]=[N+1:3]=[N-:4]",
    "[S+X3:1]-[O-:2]>>[S+0:1]=[O+0:2]",
]
for i, rxn in enumerate(REACTIONS):
    REACTIONS[i] = ReactionFromSmarts(rxn)

# files
in_file = DATA / "chembl_33.sdf.gz"
out_file = DATA / "chembl_processed.smi.gz"
num_entries = int(open(DATA / ".fetched_count").read())

# keep properties when pickling molecules
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.MolProps)


def prepare_input(mol):
    if not mol:
        return
    name = mol.GetProp("chembl_id")
    # only keep largest fragment
    mol = max(Chem.GetMolFrags(mol, asMols=True), key=lambda m: m.GetNumHeavyAtoms())
    # filter by heavy atom count
    count = mol.GetNumHeavyAtoms()
    if count > max_heavy_atoms or count < min_heavy_atoms:
        return
    # discard radicals
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return
    # sanitize
    for rxn in REACTIONS:
        apply_reaction(rxn, mol)
    err = Chem.SanitizeMol(mol, catchErrors=True)
    if err:
        return
    smi = Chem.MolToSmiles(mol)
    # roundtrip mol to smiles to mol
    roundtrip = Chem.MolFromSmiles(smi)
    if not roundtrip:
        return
    inchikey = Chem.MolToInchiKey(mol)
    return f"{smi} {name} {inchikey}\n"


count = 0
with mp.Pool(N_WORKERS) as pool, gzip.open(in_file, "r") as fi, gzip.open(
    out_file, "wt"
) as fo:
    suppl = Chem.ForwardSDMolSupplier(fi)
    for result in tqdm(
        pool.imap_unordered(prepare_input, suppl),
        total=num_entries,
        desc="Processing molecules",
    ):
        if result:
            fo.write(result)
            count += 1

print(f"Wrote '{out_file.relative_to(ROOT)}' with {count:,} entries")
(DATA / ".processed_count").write_text(str(count))
