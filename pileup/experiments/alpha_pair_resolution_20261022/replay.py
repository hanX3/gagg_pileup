from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_RESOLUTION_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'replay-source'))
if sys.argv[1] == 'fit':
    sys.argv.pop(1)
    from pileup.analysis.alpha_pair_resolution_fit import main
else:
    from pileup.analysis.alpha_pair_resolution_data import main
main()
