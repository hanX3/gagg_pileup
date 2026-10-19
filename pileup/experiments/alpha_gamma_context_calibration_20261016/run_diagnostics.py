from pathlib import Path
import sys,runpy
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'postprocessing-source/pileup/analysis'))
runpy.run_module('pileup.analysis.alpha_gamma_context_diagnostics',run_name='__main__')
