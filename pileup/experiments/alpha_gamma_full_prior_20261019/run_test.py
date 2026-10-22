from pathlib import Path
import runpy
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'));sys.path.insert(0,str(folder/'dependency-source'))
runpy.run_path(str(folder/'test-source/alpha_gamma_full_prior_test.py'),run_name='__main__')
