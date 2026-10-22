"""Run the frozen scan with project roots supplied to inherited modules."""
from pathlib import Path
import os
import sys

folder = Path(__file__).resolve().parent
project = str(folder.parents[2])
for key in ('GAGG_RESOLUTION_PROJECT', 'GAGG_PAIR_PROJECT', 'GAGG_COMPONENT_PROJECT'):
    os.environ[key] = project
sys.path.insert(0, str(folder/'replay-source'))

if __name__ == '__main__':
    if sys.argv[1] == 'fit':
        sys.argv.pop(1)
        from pileup.analysis.alpha_pair_resolution_fit import main
    else:
        from pileup.analysis.alpha_pair_resolution_data import main
    main()
