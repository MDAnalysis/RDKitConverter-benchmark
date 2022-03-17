from datetime import datetime
from pathlib import Path
from rdkit import Chem

n_jobs = -2

datapath = Path(__file__).resolve().parents[1] / "data"

timestamp = datetime.now().strftime("%Y-%m-%d")

delimiter = b"***\n"

Chem.SetDefaultPickleProperties(
    Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps)


