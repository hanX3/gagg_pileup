from pathlib import Path
import sys,runpy
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
runpy.run_module('pileup.analysis.alpha_gamma_model_selection_study',run_name='__main__')
