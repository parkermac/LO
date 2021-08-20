"""
Code to print all dt_info files.

"""

from lo_tools import Lfun, zfun
from lo_tools import hycom_functions as hfun

Ldir = Lfun.Lstart()

# initial experiment list
h_list = list(hfun.hy_dict.keys())
h_list.sort()

# specify output directory
out_dir = Ldir['data'] / 'hycom' / 'dt_lists/'

for hy in h_list:
    print('Experiment = ' + hy)
    f = open(out_dir / ('dt_info_' + hy + '.txt'), 'r')
    for line in f:
        print(' - ' + line, end='')
    f.close()
    print('\n')
