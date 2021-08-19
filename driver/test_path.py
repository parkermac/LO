"""
Code to test running python from the command line.

RESULT: we neet to use the .absolute() method in order for this
to work when called from the bash command line.

"""

import sys
from pathlib import Path

pth = Path(__file__).absolute().parent.parent / 'lo_tools' / 'lo_tools'
print(pth)

if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun, zfun
Ldir = Lfun.Lstart()
for k in Ldir.keys():
    print('%s: %s' % (k, Ldir[k]))
