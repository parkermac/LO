"""
Plot the inventory of WRF files.

RESULTS:
first good day is 2012.10.07
d4 (1.4 km grid) started in 2016
d2 (12 km) and d3 (4 km) went from 3.5 days to 4 days on 2023.11.01
still getting 3 days of d4
"""

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
import matplotlib.pyplot as plt
import pandas as pd
Ldir = Lfun.Lstart()

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.options.display.width = 0 # auto-detect full display width

in_dir = Ldir['LOo'] / 'misc'
in_fn = in_dir / 'inventory_wrf.p'
df = pd.read_pickle(in_fn)

pfun.start_plot(figsize=(20,10))
plt.close('all')

df.plot(subplots=True,ylim=(0,100),title='WRF Inventory')

plt.show()
plt.savefig(in_dir / 'inventory_wrf.png')