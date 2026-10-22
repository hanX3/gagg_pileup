"""Run the frozen inference implementation on a local observed voltage NPZ."""
from pathlib import Path
import os
import sys

folder=Path(__file__).resolve().parent
os.environ['GAGG_PAIR_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
from pileup.analysis.alpha_pair_predict import main

if __name__=='__main__':main()
