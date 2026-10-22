from pathlib import Path
import runpy
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
sys.path.insert(0,str(folder/'dependency-source'))
runpy.run_module('pileup.analysis.alpha_gamma_full_prior_study',run_name='__main__')
