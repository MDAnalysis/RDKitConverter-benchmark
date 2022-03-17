from datetime import datetime
from pathlib import Path
datapath = Path(__file__).resolve().parents[1] / "data"

timestamp = datetime.now().strftime("%Y-%m-%d")
