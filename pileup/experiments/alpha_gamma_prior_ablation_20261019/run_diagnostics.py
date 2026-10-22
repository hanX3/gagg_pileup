from pathlib import Path
import runpy
folder=Path(__file__).resolve().parent
runpy.run_path(str(folder/"postprocessing-source-v2/pileup/analysis/alpha_gamma_prior_diagnostics.py"),run_name="__main__")
