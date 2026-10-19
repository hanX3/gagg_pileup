from pathlib import Path
import runpy, sys, os
sys.dont_write_bytecode=True
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/"analysis_source"))
runpy.run_module("pileup.analysis.alpha_gamma_study",run_name="__main__")
