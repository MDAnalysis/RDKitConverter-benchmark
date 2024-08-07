import os
from datetime import datetime
from multiprocessing import cpu_count
from pathlib import Path

from MDAnalysis.converters.RDKit import _infer_bo_and_charges, _standardize_patterns
from rdkit import Chem

N_WORKERS = int(os.getenv("N_WORKERS", -1))
if N_WORKERS < 0:
    N_WORKERS = cpu_count() + 1 - N_WORKERS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

timestamp = datetime.now().strftime("%Y-%m-%d")


def apply_reaction(rxn, mol):
    """Applies a unimolecular reaction inplace. Adapted from
    ``MDAnalysis.converters.RDKit._apply_reactions``.
    """
    mol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(mol)
    while rxn.RunReactantInPlace(mol):
        mol.UpdatePropertyCache(strict=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(mol)


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
    except:
        return None
    return mol


def same_molecules(mol, ref):
    """Checks if 2 molecules are the same, using their resonance structures"""
    same = mol.HasSubstructMatch(ref)
    if not same:
        same = bool(
            Chem.ResonanceMolSupplier(mol, maxStructs=50).GetSubstructMatch(ref)
        )
    return same
