"""
This focuses on property-property plots and obs-mod plots.

It is focused on a single run and year.

In the "SET CHOICES" we pull in optional command line arguments that
allow more specific plotting choices, such as using a single obs source,
or only plotting Salish Sea stations. These are best used with
-test True when plotting by hand, as opposed to using the one_step_bottle_val_plot.py
driver.

For -test True it shows the plot on the screen instead of saving to a file.

Testing on mac:
run plot_val -gtx cas7_t0_x4b -year 2017 -test True

"""
import sys
import pandas as pd
import numpy as np
import pickle
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun, zfun, zrfun
import obsmod_functions as omfun

# command line arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas7_t1_x11ab
parser.add_argument('-otype', type=str, default='bottle') # observation type, e.g. ctd, bottle, etc.
parser.add_argument('-year', type=int) # e.g. 2019
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# more optional arguments
parser.add_argument('-single_source', type=str, default='all') # 'all' or a souce name like 'kc'
parser.add_argument('-dividing_depth', type=int, default=10) # depth [m] to delineate shallow
parser.add_argument('-do_arag', default=True, type=Lfun.boolean_string) # plot arag and pCO2
parser.add_argument('-coast', default=False, type=Lfun.boolean_string) # only plot coastal stations
parser.add_argument('-salish', default=False, type=Lfun.boolean_string) # only plot Salish stations
args = parser.parse_args()

Ldir = Lfun.Lstart()

if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

in_dir = Ldir['LOo'] / 'obsmod'

# run info
year = str(args.year)
gtx = args.gtagex
otype = args.otype

# specify input
in_fn = in_dir / ('combined_' + otype + '_' + year + '_' + gtx + '.p')
df0_dict = pickle.load(open(in_fn, 'rb'))

# add DIN field
for gtxo in df0_dict.keys():
    df0_dict[gtxo]['DIN (uM)'] = df0_dict[gtxo]['NO3 (uM)'] + df0_dict[gtxo]['NH4 (uM)']

# Get a list of available obs sources
source_list = list(df0_dict['obs']['source'].unique())

# where to put output figures
out_dir = Ldir['LOo'] / 'obsmod_val_plots'
Lfun.make_dir(out_dir)

# ========= SET CHOICES =============================================

# (1) Data choices

# Use only a single source (set to 'all' to use all available)
single_source = args.single_source # e.g. 'kc'
if single_source == 'all':
    pass
else:
    if single_source not in source_list:
        print('single source not in source_list')
        sys.exit()

# (2) Plotting choices

small = False # True for laptop size plot

H = args.dividing_depth # dividing depth [m] for deep and shallow

# Calculate  and plot values of aragonite saturation state and pCO2
do_arag = args.do_arag
# Otherwise it plots DIN and Chl

# (3) Filtering choices

fil_dict = dict() # dict to hold selected filter choices

# Set coast to True to only plot coastal stations
fil_dict['coast'] = args.coast
# Set salish to True to only plot Salish Sea stations
fil_dict['salish'] = args.salish
if fil_dict['coast'] and fil_dict['salish']:
    print('Error: Too many spatial masks!')
    sys.exit()
    
# ======== APPLY CHOICES ==================================
        
# only plot coastal stations
if fil_dict['coast']:
    for gtxo in df0_dict.keys():
        a = df0_dict[gtxo].copy()
        mask1 = (a.lat>=46) & (a.lat<49) & (a.lon>-124)
        mask2 = (a.lat>=49) & (a.lat<51) & (a.lon>-125)
        a = a.loc[(~mask1) & (~mask2),:]
        df0_dict[gtxo] = a
        
# only plot Salish Sea stations
if fil_dict['salish']:
    for gtxo in df0_dict.keys():
        a = df0_dict[gtxo].copy()
        mask1 = (a.lat>=47) & (a.lat<49) & (a.lon>-124)
        mask2 = (a.lat>=49) & (a.lat<51) & (a.lon>-125)
        mask = (~mask1) & (~mask2)
        a = a.loc[~mask,:]
        df0_dict[gtxo] = a

# start assembling some text for the plot that will include info about the filters
f_str = otype + ' ' + year + '\n' + gtx + '\n' # a string to put for info on the map
ff_str = otype + '_' + year + '_' + gtx # a string for the output .png file name

# limit which sources to use
if single_source == 'all':
    # use df_dict as-is
    f_str += 'Source = all\n'
    ff_str += '_all'
else:
    # use just one source
    f_str += 'Source = ' + single_source + '\n'
    ff_str += '_' + single_source
    for gtxo in df0_dict.keys():
        df0_dict[gtxo] = df0_dict[gtxo].loc[df0_dict[gtxo].source==single_source,:]
        
f_str += 'Dividing depth = %d m\n' % (H)

for fil in fil_dict.keys():
    if fil_dict[fil]:
        f_str += 'Filter: %s\n' % (fil)

# PLOTTING

plt.close('all')

if not do_arag:
    vn_list = ['SA','CT','DO (uM)','NO3 (uM)','NH4 (uM)','DIN (uM)',
        'DIC (uM)', 'TA (uM)', 'Chl (mg m-3)']
else:
    vn_list = ['SA','CT','DO (mg L-1)','NO3 (uM)','NH4 (uM)','pCO2 (uatm)',
        'DIC (uM)', 'TA (uM)', 'Omega']

jj_list = [1,2,3,5,6,7,9,10,11,12] # indices for the data plots

lim_dict = {'SA':(14,36),'CT':(0,20),'DO (uM)':(0,500),'DO (mg L-1)':(0,15),
    'NO3 (uM)':(0,50),'NH4 (uM)':(0,10),'DIN (uM)':(0,50),
    'DIC (uM)':(1500,2600),'TA (uM)':(1500,2600),'Chl (mg m-3)':(0,20),'Omega':(0,3),
    'pCO2 (uatm)':(0,6000)}

# create DO (mg L-1) and pCO2 (uatm)
for og in ['obs',gtx]:
    df0_dict[og]['DO (mg L-1)'] = (32 / 1000) * df0_dict[og]['DO (uM)']
if do_arag:
    import gsw
    #from PyCO2SYS import CO2SYS
    # calculate and add Aragonite Saturation State
    for og in ['obs',gtx]:
        # load data to vectors
        z = df0_dict[og]['z'].to_numpy()
        lon = df0_dict[og]['lon'].to_numpy()
        lat = df0_dict[og]['lat'].to_numpy()
        SA = df0_dict[og]['SA'].to_numpy()
        CT = df0_dict[og]['CT'].to_numpy()
        DIC0 = df0_dict[og]['DIC (uM)'].to_numpy()
        TA0 = df0_dict[og]['TA (uM)'].to_numpy()
        # Calculate derived quantities
        p = gsw.p_from_z(z, lat)
        SP = gsw.SP_from_SA(SA, p, lon, lat)
        rho = gsw.rho(SA, CT, p) # in situ density
        temp = gsw.t_from_CT(SA, CT, p) # in situ temperature
        # convert from umol/L to umol/kg using in situ dentity
        TA = 1000 * TA0 / rho
        TA[TA < 100] = np.nan
        TIC = 1000 * DIC0 / rho
        TIC[TIC < 100] = np.nan
        # See LPM/co2sys_test/test0.py for info.
        import PyCO2SYS as pyco2
        CO2dict = pyco2.sys(par1=TA, par2=TIC, par1_type=1, par2_type=2,
            salinity=SP, temperature=temp, pressure=p,
            total_silicate=50, total_phosphate=2,
            opt_pH_scale=1, opt_k_carbonic=10, opt_k_bisulfate=1)
        df0_dict[og]['Omega'] = CO2dict['saturation_aragonite']
        # also get pCO2
        df0_dict[og]['pCO2 (uatm)'] = CO2dict['pCO2']

if small:
    fs = 10
    pfun.start_plot(figsize=(13,8), fs=fs)
else:
    fs = 14
    pfun.start_plot(figsize=(20,12), fs=fs)

depth_list = ['deep', 'shallow']
c_dict = dict(zip(depth_list,['b','r']))
t_dict = dict(zip(depth_list,[.05,.15])) # vertical position of stats text

alpha = 0.3
fig = plt.figure()

ax_dict = {}
for depth_range in depth_list:
    
    df_dict = df0_dict.copy()
    
    # depth range
    if depth_range == 'shallow':
        # shallow water
        zz = -H
        for gtxo in df_dict.keys():
            df_dict[gtxo] = df_dict[gtxo].loc[df_dict[gtxo].z >= zz,:]
    elif depth_range == 'deep':
        # deep water
        zz = -H
        for gtxo in df_dict.keys():
            df_dict[gtxo] = df_dict[gtxo].loc[df_dict[gtxo].z <= zz,:]

    for ii in range(len(vn_list)):
        jj = jj_list[ii]
        if depth_range == depth_list[0]:
            ax = fig.add_subplot(3,4,jj)
            ax_dict[ii] = ax
        else:
            ax = ax_dict[ii]
        vn = vn_list[ii]           

        x = df_dict['obs'][vn].to_numpy()
        y = df_dict[gtx][vn].to_numpy()
        ax.plot(x,y,marker='.',ls='',color=c_dict[depth_range], alpha=alpha)

        if (not np.isnan(x).all()) and (not np.isnan(y).all()) and (len(x) > 0) and (len(y) > 0):
            bias = np.nanmean(y-x)
            rmse = np.sqrt(np.nanmean((y-x)**2))
            ax.text(.95,t_dict[depth_range],'bias=%0.1f, rmse=%0.1f' % (bias,rmse),
                c=c_dict[depth_range],
                transform=ax.transAxes, ha='right', fontweight='bold', bbox=pfun.bbox,
                fontsize=fs*.7,style='italic')
                            

        if jj in [9,10,11]:
            ax.set_xlabel('Observed')
        if jj in [1,5,9]:
            ax.set_ylabel('Modeled')

        # add labels to identify the model runs with the colors
        if jj == 1:
            yy = 0
            for Depth_range in c_dict.keys():
                ax.text(.05, .7 + 0.1*yy, Depth_range, c=c_dict[Depth_range],
                    transform=ax.transAxes,
                    fontweight='bold', ha='left')
                yy += 1

        ax.text(.05,.9,vn,transform=ax.transAxes, fontweight='bold')
            
        ax.axis([lim_dict[vn][0], lim_dict[vn][1], lim_dict[vn][0], lim_dict[vn][1]])
        ax.plot([lim_dict[vn][0], lim_dict[vn][1]], [lim_dict[vn][0], lim_dict[vn][1]],'-g')
        ax.grid(True)
        
# station map
ax = fig.add_subplot(1,4,4)
df_dict['obs'].plot(x='lon',y='lat',style='.g',legend=False, ax=ax)
pfun.add_coast(ax)
ax.axis([-130,-122,42,52])
pfun.dar(ax)
ax.set_xlabel('')
ax.set_ylabel('')
ax.text(.05,0,f_str,va='bottom',transform=ax.transAxes,fontweight='bold')

fig.tight_layout()

print('Plotting ' + ff_str)
sys.stdout.flush()

if args.testing:
    plt.show()
else:
    print('Saving plot to:\n %s' % (str(out_dir / (ff_str + '.png'))))
    plt.savefig(out_dir / (ff_str + '.png'))
    
