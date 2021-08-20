"""
These are snippets of python code that I often use in the LO system.
They are here for easy cut-and-paste.
"""

# This is the old way of adding the path to shared modules, which is only needed
# when operating outside of the loenv environment as we might do on mox or klone.
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent / 'lo_tools' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# This is the new way of importing the shared modules, making use of the
# "package" lo_tools
from lo_tools import Lfun, zrfun, zfun

# This gets the path to a ROMS history file on my mac
Ldir = Lfun.Lstart()
fn = (Ldir['parent'] / 'LiveOcean_roms' / 'output' /
    'cas6_v3_lo8b' / 'f2019.07.04' / 'ocean_his_0001.nc')

# Plotting stuff
import matplotlib.pyplot as plt
from cmocean import cm
from lo_tools import plotting_functions as pfun

if True:
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