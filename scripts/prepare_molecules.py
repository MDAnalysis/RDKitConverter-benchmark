import gzip
from joblib import Parallel, delayed
from tqdm.auto import tqdm
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from utils import datapath, timestamp, delimiter, n_jobs


# sanitization reactions
REACTIONS = [
    "[#6:1]-[N+:2]#[N:3]=[N-:4]>>[#6:1]-[N+0:2]=[N+1:3]=[N-:4]",
    "[S+X3:1]-[O-:2]>>[S+0:1]=[O+0:2]",
]
for i, rxn in enumerate(REACTIONS):
    REACTIONS[i] = ReactionFromSmarts(rxn)

in_file = next(datapath.glob("chembl_*_filtered.smi.gz"))
out_file = datapath / f"chembl_{timestamp}_prepared.rdkit.gz"
num_entries = int(open(datapath / f".chembl_filtered_count").read())


def prepare_input(smi):
    mol = Chem.MolFromSmiles(smi[:-1], sanitize=False)
    name = mol.GetProp("_Name")
    # only keep largest fragment
    mol = max(Chem.GetMolFrags(mol, asMols=True), 
              key=lambda m: m.GetNumHeavyAtoms())
    # sanitize
    for rxn in REACTIONS:
        while rxn.RunReactantInPlace(mol):
            pass
    err = Chem.SanitizeMol(mol, catchErrors=True)
    if err:
        return
    mol.SetProp("_Name", name)
    bytestr = mol.ToBinary() + delimiter
    return bytestr

with gzip.open(in_file, "rt") as smiles_fi, \
     gzip.open(out_file, "wb") as fo:

    count = 0
    for bytestr in Parallel(n_jobs=n_jobs)(delayed(prepare_input)(smi)
                            for smi in tqdm(smiles_fi, total=num_entries)):
        if bytestr:
            fo.write(bytestr)
            count += 1
(datapath / ".molecules_count").write_text(str(count))
