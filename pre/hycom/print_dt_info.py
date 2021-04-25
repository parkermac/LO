"""
Code to print all dt_info files.

"""

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import zfun
Ldir = Lfun.Lstart()

import pickle

import hfun
from importlib import reload
reload(hfun)

# initial experiment list
h_list = list(hfun.hy_dict.keys())
h_list.sort()

# specify output directory
out_dir = Ldir['data'] / 'hycom' / 'dt_lists/'

for hy in h_list:
    print('Experiment = ' + hy)
    f = open(out_dir / ('dt_info_' + hy + '.txt'), 'r')
    for line in f:
        sys.stdout.write(' - ' + line)
    f.close()
    print('\n')
