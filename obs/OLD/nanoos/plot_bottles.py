"""
Code to plot bottle data in informative ways.
"""

from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd

import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# Choices
source = 'nanoos'
otype = 'bottle'
year = 2021

# file names
out_dir = Ldir['LOo'] / 'obs' / source / otype
ys = str(year)
out_fn = out_dir / (str(year) + '.p')
info_out_fn = out_dir / ('info_' + str(year) + '.p')

# load data
df = pd.read_pickle(out_fn)
info_df = pd.read_pickle(info_out_fn)

# create in situ density
SA = df.SA.to_numpy()
CT = df.CT.to_numpy()
lon = df.lon.to_numpy()
lat = df.lat.to_numpy()
z = df.z.to_numpy()
p = gsw.p_from_z(z, lat)
rho = gsw.rho(SA, CT, p) # [kg m-3] not anomaly

# Variables
DO_uM = df['DO (uM)'].to_numpy()
NO3_uM = df['NO3 (uM)'].to_numpy()
NH4_uM = df['NH4 (uM)'].to_numpy()

# Calculate what DO would be if it was saturated
# DO_umol_kg = df['DO (uM)'] * 1000 / rho
DOsol_umol_kg = gsw.O2sol(SA, CT, p, lon, lat)
DOsol_uM = DOsol_umol_kg * rho / 1000
AOU = DOsol_uM - DO_uM # Apparent Oxygen Utilization

# Calculate what DO would be for our boundary conditions
Ofun_bio = Lfun.module_from_file('Ofun_bio', Ldir['LO'] / 'forcing' / 'ocn00' / 'Ofun_bio.py')
SP = gsw.SP_from_SA(SA, p, lon, lat)
DObc_uM = Ofun_bio.create_bio_var(SP, 'oxygen')

# plotting
plt.close('all')
pfun.start_plot(figsize=(12,7))

m_dict = {'apr':(df.cruise == 'RC0051').to_numpy(),
        'jul':(df.cruise == 'RC0058').to_numpy(),
        'sept':(df.cruise == 'RC0063').to_numpy()}

c_dict = {'apr':'g', 'jul':'gold', 'sept':'r'}

fig, axes = plt.subplots(nrows=2, ncols=2)

for m in m_dict.keys():
    
    ax = axes[0,0]
    ax.plot(DO_uM[m_dict[m]], z[m_dict[m]], linestyle='', marker='.', color=c_dict[m])
    ax.grid(True)
    ax.set_title('DO (uM)')
    
    ax = axes[0,1]
    ax.plot(AOU[m_dict[m]], z[m_dict[m]], linestyle='', marker='.', color=c_dict[m])
    ax.grid(True)
    ax.set_title('AOU (uM)')
    
    ax = axes[1,0]
    ax.plot(NO3_uM[m_dict[m]], z[m_dict[m]], linestyle='', marker='.', color=c_dict[m])
    ax.grid(True)
    ax.set_title('NO3 (uM)')
    
    ax = axes[1,1,]
    ax.plot(NH4_uM[m_dict[m]], z[m_dict[m]], linestyle='', marker='.', color=c_dict[m])
    ax.grid(True)
    ax.set_title('NH4 (uM)')

plt.show()


