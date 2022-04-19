"""
Code to test that loenv works, in particular that the lo_tools package
is present and functioning.

"""

from lo_tools import Lfun

Ldir = Lfun.Lstart()
for k in Ldir.keys():
    print('%s: %s' % (k, Ldir[k]))
