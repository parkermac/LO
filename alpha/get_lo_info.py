"""
This is the one place where you set the path structure of the LO code.
The info is stored in the dict Ldir.

All paths are stored as pathlib.Path objects

This program is meant to be loaded as a module by Lfun which then adds more
entries to the Ldir dict based on which model run you are workgin on.

Users should not edit this directly but instead should copy it to
LO_user/alpha/user_get_lo_info.py where LO_user is at the same level as LO.
Then they can edit that version, and put LO_user in their own GitHub repo.
"""
import os
from pathlib import Path

# defaults that should work on all machines
parent = Path(__file__).absolute().parent.parent.parent
LO = parent / 'LO'
LOo = parent / 'LO_output'
LOu = parent / 'LO_user'
data = parent / 'LO_data'
roms = parent / 'LO_roms'
roms1 = parent / 'BLANK'
roms2 = parent / 'BLANK'

# default for linux machines
which_matlab = '/usr/local/bin/matlab'

HOME = Path.home()
try:
    HOSTNAME = os.environ['HOSTNAME']
except KeyError:
    HOSTNAME = 'BLANK'

if str(HOME) == '/Users/pm8':
    lo_env = 'pm_mac'
    roms = parent / 'LiveOcean_roms'
    which_matlab = '/Applications/MATLAB_R2020a.app/bin/matlab'

elif (str(HOME) == '/home/parker') & ('boiler' in HOSTNAME):
    lo_env = 'pm_boiler'
    roms1 = Path('/pgdat1/parker/LO_roms')

elif (str(HOME) == '/home/parker') & ('perigee' in HOSTNAME):
    lo_env = 'pm_perigee'
    roms1 = Path('/boildat1/parker/LO_roms')
    roms2 = Path('/data2/parker/LO_roms')
  
Ldir = dict()
Ldir['lo_env'] = lo_env
Ldir['parent'] = parent
Ldir['LO'] = LO
Ldir['LOo'] = LOo
Ldir['LOu'] = LOu
Ldir['data'] = data
Ldir['roms'] = roms
Ldir['roms1'] = roms1
Ldir['roms2'] = roms2
Ldir['which_matlab'] = which_matlab

