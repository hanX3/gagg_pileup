from pathlib import Path
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
from analysis import alpha_gamma_short_window_native as module
module.PROJECT=folder.parents[2]
module.main()
