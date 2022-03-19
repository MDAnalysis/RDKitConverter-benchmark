# RDKitConverter benchmark

This repository benchmarks the ability of MDAnalysis' `RDKitConverter` to infer bond orders and charges from molecules with all hydrogens explicit.

## Instructions

Requires conda (or mamba) on a Linux machine.

Start by cloning this repository:
```shell
git clone https://github.com/MDAnalysis/RDKitConverter-benchmark.git
```

Then install the python dependencies with `make install`:
```shell
# to speed things up you can use mamba:
make install CONDA=mamba
```
This will create a separate conda environment called `rdkitconverter`.

Finally, run the benchmark:
```shell
make
```

Run `make help` to get a list of available commands.

The results are available in the `results/` directory:
- `results.json`, a JSON file listing all the necessary information
- `failed_molecules.smi`, a SMILES file containing the molecules that failed the test
- `failed_molecules.html`, an interactive table displaying the failed molecules

## Methods

The benchmark will fetch ChEMBL 30 as an SDF file and process the molecules the following way:
- Discard molecules that could not be read or sanitized by RDKit
- Keep only the largest fragment
- Keep only molecules with 2 to 50 heavy atoms
- Discard molecules with radicals
- Drop duplicate molecules based on their InchiKey

Once the data is fetched and standardized, the benchmark can be run. The benchmark will start by preparing a "reduced" version of the molecule by removing bond orders and formal charges. This is done to mimic the minimal information available in most topology files for MD simulations.

The RDKitConverter might give different results depending on the order of atoms in the molecule. For that reason, the benchmark will enumerate reordered version of the molecule so that each atom appears in the first position once.  
This is done by reading a SMILES of the molecule rooted at the given atom, so that the other atoms of the molecule are reordered in a realistic way.

Finally, the reordered "reduced" molecule goes through the MDAnalysis code responsible for inferring bond orders and formal charges.

During the enumeration of reordered molecules, if any of the inferred molecules fails to match with the original molecule or one of its resonance structures, the whole test fails for that molecule.

Prototype code:
```python
def benchmark_mol(reference_mol):
    reduced_mol = remove_bond_orders_and_charges(reference_mol)
    for mol in enumerate_reordered_mol(reduced_mol):
        mol = infer_bond_orders_and_charges(mol)
        valid = is_same_mol_or_resonance_structure(mol, reference_mol)
        if not valid:
            return "FAILED"
    return "SUCCESS"
```
