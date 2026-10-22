from pathlib import Path
import runpy,sys
folder=Path(__file__).resolve().parent
sys.path[:0]=[str(folder/"source"),str(folder/"dependency-source")]
runpy.run_module("pileup.analysis.alpha_gamma_prior_study",run_name="__main__")
