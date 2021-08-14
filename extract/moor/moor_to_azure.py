"""
Code to simplify pushing a mooring extraction to azure.
"""
from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart()

# choose the file
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'moor'
moor_name = Lfun.choose_item(in_dir, tag='.nc', exclude_tag='', itext='** Choose mooring extraction from list **')
moor_fn = in_dir / moor_name

# this will give screen output of the URL for sharing
Lfun.copy_to_azure(moor_fn, moor_name, 'pm_share', Ldir)