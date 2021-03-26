"""
Code to test running python from the command line.

RESULT: we neet to use the .absolute() method in order for this
to work when called from the bash command line.

"""

import sys
from pathlib import Path

pth = Path(__file__).absolute().parent.parent / 'alpha'
print(pth)

if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun, zfun
