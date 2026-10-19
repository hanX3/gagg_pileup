"""Rebuild post-test figures/tables using the saved reporting source."""
from pathlib import Path
import argparse
import importlib.util
import sys

folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
sys.path.insert(1,str(folder.parents[2]))
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args()
spec=importlib.util.spec_from_file_location('pileup.analysis.frozen_postprocessing',
    folder/'postprocessing/alpha_gamma_regularized_report.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.summarize(folder,args.artifacts)
module.make_figures(folder,args.artifacts)
module.comparison_tables(args.artifacts)
