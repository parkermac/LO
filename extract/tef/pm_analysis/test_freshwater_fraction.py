"""
Code to test the sensitivity of freshwater flux to the choice of Socn.

RESULT: Qfw is weakly dependent on Socn, but I think it is best to use
Socn = tef_df['salt_in'].max() from the mouth of the system.

Socn = 33.75 in this example.

"""
from pathlib import Path
import sys
pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import flux_fun

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_Qfw(tef_df, Socn):
    # get mean freshwater flux for a given Socn
    Qfw = (-( tef_df['Qin']*(Socn-tef_df['salt_in']) +
        tef_df['Qout']*(Socn-tef_df['salt_out']) ).mean()/Socn)
    return Qfw

# specify section and bulk folder
in_dir = Path('/Users/pm8/Documents/LO_output/extract/cas6_v3_lo8b/tef/bulk_2018.01.01_2018.12.31')

# get Qfw vs. Socn over a range of Socn at ai1
sect_name = 'ai1'
gridname = 'cas6'
tef_df, _, _, _ = flux_fun.get_two_layer(in_dir, sect_name, gridname)
# get Qfw over a reange of Socn at this section
Qfw_ser = pd.Series(index = np.arange(15, 41))
for Socn in np.arange(15, 41):
    Qfw_ser[Socn] = get_Qfw(tef_df, Socn)
    
# Get Socn at the mouth of the system, jdf1
sect_name = 'jdf1'
gridname = 'cas6'
tef_df1, _, _, _ = flux_fun.get_two_layer(in_dir, sect_name, gridname)
Socn = tef_df1['salt_in'].max()
Qfw = get_Qfw(tef_df, Socn)

plt.close('all')
ax = Qfw_ser.plot()
ax.plot(Socn, Qfw, '*r')
ax.set_ylim(bottom=0)
ax.grid(True)
ax.set_xlabel('Socn')
ax.set_ylabel('Qfw [m3/s]')
plt.show()