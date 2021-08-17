"""
These are snippets of python code that I often use in the LO system.
They are here for easy cut-and-paste.
"""

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
import Lfun
import zrfun
import zfun

Ldir = Lfun.Lstart()
fn = (Ldir['parent'] / 'LiveOcean_roms' / 'output' /
    'cas6_v3_lo8b' / 'f2019.07.04' / 'ocean_his_0001.nc')

import matplotlib.pyplot as plt
import plotting_functions as pfun
from cmocean import cm

if False:
    # PLOTTING
    plt.close('all')
    pfun.start_plot(figsize=(8,12))
    fig = plt.figure()

    ax = fig.add_subplot(111)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    #ax.set_axis_off()

    fig.tight_layout()
    plt.show()
    pfun.end_plot()