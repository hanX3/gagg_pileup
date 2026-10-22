"""Use this experiment's saved Python sources and the unchanged local inputs."""
import os
from pathlib import Path
import sys

project = Path(__file__).resolve().parents[3]
os.environ['GAGG_ML_REFRESH_PROJECT'] = str(project)
sys.path.insert(0, str(Path(__file__).resolve().parent / 'inference-source'))
from pileup.analysis.alpha_ml_refresh import main

if __name__ == '__main__':
    main()
