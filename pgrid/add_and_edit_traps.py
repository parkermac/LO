"""
Code to do traps_placement and edit things by hand.
"""

from lo_tools import Lfun

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',type=str) # e.g. oly2
args = parser.parse_args()

Ldir = Lfun.Lstart(gridname=args.gridname)

# User edits
traps_pre_version = 'trapsP01'
traps_list = ['moh20_triv','moh20_wwtp','was24_wwtp'] # from LO/pre/trapsP01 2025.08.27
exclude_list = []
if Ldir['gridname'] == 'oly2':
    exclude_list = ['Lynch Cove', 'ALDERBROOK RESORT & SPA', 'Tahuya']
# end User edits

print('\nRunning traps_placement and excluding selected items')
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import pandas as pd
cmd_list = ['python', str(Ldir['LO'] / 'pre' / traps_pre_version / 'traps_placement.py'), '-g', Ldir['gridname']]
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if len(stderr) > 0:
    print(stderr.decode())

# then edit out excluded triv and wwtp
for traps in traps_list:
    traps_fn = Ldir['grid'] / (traps + '_info.csv')
    if traps_fn.is_file():
        traps_df = pd.read_csv(traps_fn, index_col='rname')
        for item in exclude_list:
            try:
                traps_df = traps_df.drop(item)
                print(' - Excluding %s from %s' % (item, traps))
            except KeyError:
                pass
        traps_df.to_csv(traps_fn)