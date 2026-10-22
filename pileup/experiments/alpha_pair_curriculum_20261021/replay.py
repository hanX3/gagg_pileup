from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_PAIR_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
from pileup.analysis.alpha_pair_study import main
main()
