"""Rebuild tables and export figures from the saved post-test reporting code."""
from pathlib import Path
import argparse
import importlib.util
import sys

FOLDER = Path(__file__).resolve().parent
sys.path.insert(0, str(FOLDER/'source'))
sys.path.insert(1, str(FOLDER.parents[2]))


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, FOLDER/'postprocessing'/filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--previous', type=Path, required=True)
args = parser.parse_args()
load('pileup.analysis.alpha_gamma_regularized_report', 'alpha_gamma_regularized_report.py')
report = load('pileup.analysis.alpha_gamma_model_selection_report', 'alpha_gamma_model_selection_report.py')
report.summarize(FOLDER, args.artifacts)
report.figures(FOLDER, args.artifacts, args.previous)
