"""
Code to test that we are getting the right environment, aimed at
making sure crontab works.

"""

from pathlib import Path
import sys

pth = Path(__file__).absolute().parent.parent / 'lo_tools' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))

import Lfun
Ldir = Lfun.Lstart()
for k in Ldir.keys():
    print('%s: %s' % (k, Ldir[k]))
