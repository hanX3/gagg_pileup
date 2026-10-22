from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(folder.parents[2])
sys.path.insert(0,str(folder/'postprocessing-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_resolution_report import main
    main()
