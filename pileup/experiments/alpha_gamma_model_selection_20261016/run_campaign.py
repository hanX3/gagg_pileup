from pathlib import Path
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
import pileup.analysis.alpha_gamma_regularized_campaign as module
module.PROJECT=folder.parents[2]
module.main()
