import os
from datetime import datetime
from pathlib import Path
from multiprocessing import cpu_count
from rdkit import Chem
from MDAnalysis.converters.RDKit import (_infer_bo_and_charges,
                                         _standardize_patterns)


N_WORKERS = int(os.getenv("N_WORKERS", -1))
if N_WORKERS < 0:
    N_WORKERS = cpu_count() + 1 - N_WORKERS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

timestamp = datetime.now().strftime("%Y-%m-%d")


def apply_reaction(rxn, mol):
    """Adapted from MDAnalysis.converters.RDKit._run_reaction"""
    for _ in range(mol.GetNumAtoms()):
        mol.UpdatePropertyCache(strict=False)
        Chem.Kekulize(mol)
        products = rxn.RunReactants((mol,))
        if products:
            # structure: tuple[tuple[mol]]
            # length of 1st tuple: number of matches in reactant
            # length of 2nd tuple: number of products yielded by the reaction
            # if there's no match in reactant, returns an empty tuple
            # here we have at least one match, and the reaction used yield
            # a single product hence `products[0][0]`
            product = products[0][0]
            # apply the next reaction to the product
            mol = product
        # exit the loop if there was nothing to transform
        else:
            break
    mol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(mol)
    return mol


def add_Hs_remove_bo_and_charges(mol):
    """Add hydrogens and remove bond orders and charges from a molecule"""
    mH = Chem.AddHs(mol)
    for atom in mH.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
    for bond in mH.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    mH.UpdatePropertyCache(strict=False)
    return mH


def enumerate_reordered_mol(mol):
    """Enumerates all possible starting atoms for a given molecule"""
    for root_atom in mol.GetAtoms():
        smi = Chem.MolToSmiles(mol, rootedAtAtom=root_atom.GetIdx())
        reordered_mol = Chem.MolFromSmiles(smi, sanitize=False)
        for atom in reordered_mol.GetAtoms():
            atom.SetNoImplicit(True)
        reordered_mol.UpdatePropertyCache(strict=False)
        yield reordered_mol


def assign_bond_orders_and_charges(mol):
    """Returns a sanitized molecule with infered bond orders and charges"""
    try:
        _infer_bo_and_charges(mol)
        mol = _standardize_patterns(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return None
    return mol


def same_molecules(mol, ref):
    """Checks if 2 molecules are the same, using their resonance structures"""
    same = mol.HasSubstructMatch(ref)
    if not same:
        same = bool(Chem.ResonanceMolSupplier(mol, maxStructs=50)
                    .GetSubstructMatch(ref))
    return same
