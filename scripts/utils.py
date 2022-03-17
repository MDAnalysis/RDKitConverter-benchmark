from datetime import datetime
from pathlib import Path
from rdkit import Chem
from MDAnalysis.converters.RDKit import _infer_bo_and_charges, _standardize_patterns

n_jobs = -2

datapath = Path(__file__).resolve().parents[1] / "data"

timestamp = datetime.now().strftime("%Y-%m-%d")

delimiter = b"***\n"

Chem.SetDefaultPickleProperties(
    Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps)


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
    # go through each possible starting atom
    for root_atom in mol.GetAtoms():
        smi = Chem.MolToSmiles(mol, rootedAtAtom=root_atom.GetIdx())
        reordered_mol = Chem.MolFromSmiles(smi, sanitize=False)
        for atom in reordered_mol.GetAtoms():
            atom.SetNoImplicit(True)
        reordered_mol.UpdatePropertyCache(strict=False)
        yield reordered_mol


def assign_bond_orders_and_charges(mol):
    """Returns a sanitized molecule with infered bond orders and charges"""
    _infer_bo_and_charges(mol)
    mol = _standardize_patterns(mol)
    Chem.SanitizeMol(mol)
    return mol


def is_isomorphic(mol, ref, useChirality=False):
    return (mol.HasSubstructMatch(ref, useChirality=useChirality)
            and ref.HasSubstructMatch(mol, useChirality=useChirality))


def is_isomorphic_or_resonance_structure(mol, ref):
    """Checks if 2 molecules are equal. If not, enumerates the resonance
    structures of the first molecule and checks again"""
    isomorphic = is_isomorphic(mol, ref)
    if not is_isomorphic:
        # try resonance structures
        isomorphic = any(is_isomorphic(res_mol, ref)
                         for res_mol in Chem.ResonanceMolSupplier(mol))
    return isomorphic
