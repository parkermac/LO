"""
Code to read the ROMS text file varinfo.dat to construct the variable
dimensions and attributes used when writing to NetCDF.

This works but is an incomplete solution.  For now it is just an experiment.
"""

from lo_tools import Lfun

Ldir = Lfun.Lstart()

fn = Ldir['parent'] / 'LiveOcean_roms' / 'LO_ROMS' / 'ROMS' / 'External' / 'varinfo.dat'

vn_list = ['salt', 'temp', 'u', 'v', 'zeta', 'ubar', 'vbar']

flag = 0
with open(fn, 'r') as f:
    for a in f:
        
        if flag==1:
            flag = 0
            try:
                long_name = a.split("'")[1]
                if 'climatology' in long_name:
                    print('%s: %s' % (vn, long_name))
                flag = 0
            except IndexError:
                pass
        if a[0] == "'" and flag==0:
            vn = a.split("'")[1]
            if vn in vn_list:
                flag = 1
                
