from pathlib import Path
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder/'inference-source'))
from analysis import alpha_gamma_shared_tau_native as module
module.PROJECT = folder.parents[2]
if __name__ == '__main__':
    module.main()
