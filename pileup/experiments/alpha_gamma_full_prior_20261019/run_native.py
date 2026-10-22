from pathlib import Path
import argparse
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'));sys.path.insert(0,str(folder/'dependency-source'))
from pileup.analysis import alpha_gamma_regularized_campaign as campaign
p=argparse.ArgumentParser()
p.add_argument('--stage',choices=['preflight','test_mixed'],required=True)
p.add_argument('--inputs',type=Path,required=True)
p.add_argument('--workers',type=int,default=4)
a=p.parse_args();campaign.PROJECT=folder.parents[2]
campaign.run(folder,a.stage,a.inputs,a.workers)
