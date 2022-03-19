import gzip
import multiprocessing as mp
from tqdm.auto import tqdm
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from utils import ROOT, DATA, N_WORKERS, apply_reaction


# sanitization reactions
REACTIONS = [
    "[#6:1]-[N+:2]#[N:3]=[N-:4]>>[#6:1]-[N+0:2]=[N+1:3]=[N-:4]",
    "[S+X3:1]-[O-:2]>>[S+0:1]=[O+0:2]",
]
for i, rxn in enumerate(REACTIONS):
    REACTIONS[i] = ReactionFromSmarts(rxn)

in_file = DATA / "chembl_fetched.smi.gz"
out_file = DATA / "chembl_processed.smi.gz"
num_entries = int(open(DATA / ".fetched_count").read())


def prepare_input(line):
    smi, name = line[:-1].split(" ")
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    # only keep largest fragment
    mol = max(Chem.GetMolFrags(mol, asMols=True), 
              key=lambda m: m.GetNumHeavyAtoms())
    # discard radicals
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return
    # sanitize
    for rxn in REACTIONS:
        mol = apply_reaction(rxn, mol)
    err = Chem.SanitizeMol(mol, catchErrors=True)
    if err:
        return
    smi = Chem.MolToSmiles(mol)
    inchikey = Chem.MolToInchiKey(mol)
    return f"{smi} {name} {inchikey}\n"


count = 0
with mp.Pool(N_WORKERS) as pool, \
     gzip.open(in_file, "rt") as fi, \
     gzip.open(out_file, "wt") as fo:

    for result in tqdm(pool.imap_unordered(prepare_input, fi),
                       total=num_entries, desc="Cleaning molecules"):
        if result:
            fo.write(result)
            count += 1

print(f"Wrote '{out_file.relative_to(ROOT)}' with {count:,} entries")
(DATA / ".processed_count").write_text(str(count))
