from pathlib import Path
import runpy,sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
sys.path.insert(1,str(folder.parents[2]))
runpy.run_module('pileup.analysis.alpha_gamma_regularized_study',run_name='__main__')
