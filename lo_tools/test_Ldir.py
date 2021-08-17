"""
Code to test getting Ldir using the new lo_tools module.

"""

from lo_tools import Lfun

Ldir = Lfun.Lstart()

print(' Contents of Ldir '.center(60,'-'))
for k in Ldir.keys():
    print('%s: %s' % (k, str(Ldir[k])))
