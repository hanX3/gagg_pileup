from pathlib import Path
import argparse
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder / 'inference-source'))
from analysis import alpha_gamma_advanced_native_local as module
module.PROJECT = folder.parents[2]
p = argparse.ArgumentParser()
p.add_argument('--workers', type=int, default=4)
a = p.parse_args()
module.run(folder, workers=a.workers)
