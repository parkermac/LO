"""
Script that maps TRAPS to grid.
Takes in lat/lon coordinates of TRAPS (generated by the ecology_excel2netcdf.py file)
Then creates triv_info.csv and wwtp_info.csv with grid indices (to place sources)
and river direction.

Option to plot source location for review as well

To run (from ipython):
run traps_placement.py -g cas7

To look at plots, run:
run traps_placement.py -g cas7 -test True

To change the Ecology data version, run from ipython with:
run traps_placement.py -tD trapsD##
(default is trapsD00)
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

import xarray as xr
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import cmocean
import argparse
import sys
from lo_tools import plotting_functions as pfun
import traps_helper
from lo_tools import Lfun

#################################################################################
#                              Argument parsing                                 #
#################################################################################

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str) 
# -test True will output plots
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# add ctag
parser.add_argument('-ctag', type=str, default='lo_base')
# add Ecology data version
parser.add_argument('-tD', type=str, default='trapsD00')
args = parser.parse_args()
ctag = args.ctag
trapsD = args.tD

# make sure grid name was provided
argsd = args.__dict__
if argsd['gridname'] == None:
    print('*** Missing required argument to forcing_argfun.intro(): gridname')
    sys.exit()

Ldir = Lfun.Lstart(gridname=args.gridname)

#################################################################################
#                     Function to place TRAPS on grid                           #
#################################################################################

def traps_placement(source_type,Ldir):
    '''
    Function that looks at all of the rivers and marine point sources 
    in SSM, and identified where to place these tiny rivers and point sources (TRAPS)
    in the LiveOcean grid. 

    Input:
        source_type:    'riv' or 'wwtp'

    Output:
        wwtp_info.csv or triv_info.csv saved in LO_data/grids/[gridname]
        Figures with source locations saved in LO_data/grids/[gridname]
    '''

    # open data if dealing with either river or wwtp
    if source_type == 'wwtp':
        output_fn = 'wwtp_info.csv'
        # read data
        wwtp_fn = Ldir['data'] / trapsD / 'all_point_source_data.nc'
        source_ds = xr.open_dataset(wwtp_fn)
    elif source_type == 'riv':
        output_fn = 'triv_info.csv'
        # read data
        riv_fn = Ldir['data'] / trapsD / 'all_nonpoint_source_data.nc'
        source_ds = xr.open_dataset(riv_fn)

    # get the grid data
    grid_fn = Ldir['grid'] / 'grid.nc'
    print(grid_fn)
    ds = xr.open_dataset(grid_fn)
    z = -ds.h.values
    mask_rho = np.transpose(ds.mask_rho.values)
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    X = lon[0,:] # grid cell X values
    Y = lat[:,0] # grid cell Y values

    # read overlapping rivers
    repeatrivs_fn = Ldir['data'] / trapsD / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)

    # initialize dataframe to save results
    rowcol_df = pd.DataFrame()
    SSMrivll_df = pd.DataFrame()
    snames = []
    ringnums = []

    # get list of source names
    snames = source_ds.name.values

    # loop through all of the traps
    for source in snames:

        # get lat/lon of source
        x = source_ds.lon[source_ds.name == source].values[0]
        y = source_ds.lat[source_ds.name == source].values[0]
        
        # place tiny rivers
        if source_type == 'riv':

            # check if river already in LiveOcean
            SSM_repeats = repeatrivs_df['SSM_rname'] # get names of repeat rivers
            # don't add rivers that already exist
            # Also drop Willamette R
            if (SSM_repeats.str.contains(source).any()) or (source == 'Willamette R'):
                continue 
            else:
                # add river to LiveOcean if not pre-existing
                # get nearest coastal gridcell
                [row,col,idir,isign,uv,ringnum] = traps_helper.get_nearest_coastal_cell_riv(source,x,y,X,Y,mask_rho)
                # save coordinates to rowcol_df
                rowcol_df.loc[source, 'row_py'] = row 
                rowcol_df.loc[source, 'col_py'] = col
                rowcol_df.loc[source, 'idir']   = idir
                rowcol_df.loc[source, 'isign']  = isign
                rowcol_df.loc[source, 'uv']     = uv
                SSMrivll_df.loc[source,'Lon'] = x
                SSMrivll_df.loc[source,'Lat'] = y
                # save ringnum for tracking. ringnum = zero means the source didn't move
                ringnums = ringnums + [ringnum]
        
        # place point sources
        else: # wwtps are easier
            # get nearest coastal gridcell
            [row,col,ringnum] = traps_helper.get_nearest_coastal_cell_wwtp(source,x,y,X,Y,mask_rho)
            # save coordinates to rowcol_df
            rowcol_df.loc[source, 'row_py'] = row 
            rowcol_df.loc[source, 'col_py'] = col
            # save ringnum for tracking. ringnum = zero means the source didn't move
            ringnums = ringnums + [ringnum]

        # update name of index column
        rowcol_df.index.name = 'rname'

#################################################################################
#                                 Initialize plot                               #
#################################################################################

    if args.testing == True:

        plon, plat = pfun.get_plon_plat(lon,lat)

        # make a version of z with nans where masked
        zm = z.copy()
        zm[np.transpose(mask_rho) == 0] = np.nan
        zm[np.transpose(mask_rho) != 0] = -1

        # add water and grid lines
        fig = plt.figure(figsize=(5,9))
        ax = fig.add_subplot(111)
        pfun.add_coast(ax,color='black')
        ax.pcolormesh(plon, plat, zm, edgecolor='aliceblue', linewidth=0.5,
                      vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        # Full grid
        ax.set_xlim(-130,-121.5)
        ax.set_ylim(42,52)

        # format plot
        pfun.dar(ax) 
        ax.set_title('placement of {}s'.format(source_type), fontsize = 18)
        ax.set_ylabel('Lat', fontsize = 16)
        ax.set_xlabel('Lon', fontsize = 16)
        ax.tick_params(axis='both', which='major', labelsize=14)
        plt.locator_params(nbins=4)

#################################################################################
#                  Plot tiny rivers and pre-existing LO rivers                  #
#################################################################################

        if source_type == 'riv':
            lon_u = ds.lon_u.values
            lat_u = ds.lat_u.values
            lon_v = ds.lon_v.values
            lat_v = ds.lat_v.values

            # label counters
            first_label_triv = True
            first_label_SSM = True
            first_label_LO = True

            # get pre-exiting LO river information
            LOrivs_fn = Ldir['grid'] / 'river_info.csv'
            LOrivs_df = pd.read_csv(LOrivs_fn)

            # Tiny river plotting
            for i,rn in enumerate(rowcol_df.index):
                # These are indices (python, zero-based) into either the
                # u or v grids.
                ii = rowcol_df.loc[rn,'col_py']
                jj = rowcol_df.loc[rn,'row_py']
                
                uv = rowcol_df.loc[rn,'uv']
                isign = rowcol_df.loc[rn,'isign']
                idir = rowcol_df.loc[rn,'idir']

                # plot original sources
                if first_label_SSM:
                    ax.plot(SSMrivll_df['Lon'][i],SSMrivll_df['Lat'][i], linestyle='none', markeredgecolor ='k',
                        marker='D', color='yellowgreen', label='Original lat/lon coordinates')
                    first_label_SSM = False
                else:
                    ax.plot(SSMrivll_df['Lon'][i],SSMrivll_df['Lat'][i], linestyle='none', markeredgecolor ='k',
                        marker='D', color='yellowgreen')
                    
                # remove rivers that are outside of the LO grid
                if np.isnan(ii):
                    continue
                else:
                    ii = int(ii)
                    jj = int(jj)

                # plot the new sources in LO grid (accounting for river shift)
                if uv == 'u' and isign == 1:
                    # River source on W side of rho cell
                    ax.plot(lon_u[jj,ii], lat_u[jj,ii],'>r')
                    if first_label_triv:
                        ax.plot(lon[jj,ii+1], lat[jj,ii+1],color='purple', marker='o',
                         linestyle = 'None', label='Placed location in LO')
                        first_label_triv = False
                    else:
                        ax.plot(lon[jj,ii+1], lat[jj,ii+1],color='purple', marker='o', linestyle = 'None')
                    # add river name
                    ax.text(lon[jj,ii+1], lat[jj,ii+1]+0.003, rn, color = 'purple',
                            fontsize=12, horizontalalignment='center')
                    # plot tracks
                    ax.plot([SSMrivll_df['Lon'][i], lon[jj,ii+1]],
                            [SSMrivll_df['Lat'][i], lat[jj,ii+1]],
                            color='hotpink', linewidth=0.5)
                elif uv == 'u' and isign == -1:
                    # River source on E side of rho cell
                    ax.plot(lon_u[jj,ii], lat_u[jj,ii],'<r')
                    ax.plot(lon[jj,ii], lat[jj,ii],color='purple', marker='o', linestyle = 'None')
                    # add river name
                    ax.text(lon[jj,ii], lat[jj,ii]+0.003, rn, color = 'purple',
                            fontsize=12, horizontalalignment='center')
                    # plot tracks
                    ax.plot([SSMrivll_df['Lon'][i], lon[jj,ii]],
                            [SSMrivll_df['Lat'][i], lat[jj,ii]],
                            color='hotpink', linewidth=0.5)
                elif uv == 'v' and isign == 1:
                    # River source on S side of rho cell
                    ax.plot(lon_v[jj,ii], lat_v[jj,ii],'^b')
                    ax.plot(lon[jj+1,ii], lat[jj+1,ii],color='purple', marker='o', linestyle = 'None')
                    # add river name
                    ax.text(lon[jj+1,ii], lat[jj+1,ii]+0.003, rn, color = 'purple',
                            fontsize=12, horizontalalignment='center')
                    # plot tracks
                    ax.plot([SSMrivll_df['Lon'][i], lon[jj+1,ii]],
                            [SSMrivll_df['Lat'][i], lat[jj+1,ii]],
                            color='hotpink', linewidth=0.5)
                elif uv == 'v' and isign == -1:
                    # River source on N side of rho cell
                    ax.plot(lon_v[jj,ii], lat_v[jj,ii],'vb')
                    ax.plot(lon[jj,ii], lat[jj,ii],color='purple', marker='o', linestyle = 'None')
                    # add river name
                    ax.text(lon[jj,ii], lat[jj,ii]+0.003, rn, color = 'purple',
                            fontsize=12, horizontalalignment='center')
                    # plot tracks
                    ax.plot([SSMrivll_df['Lon'][i], lon[jj,ii]],
                            [SSMrivll_df['Lat'][i], lat[jj,ii]],
                            color='hotpink', linewidth=0.5)

            # Pre-existing LO river plotting
            for rn in LOrivs_df.index:
                # These are indices (python, zero-based) into either the
                # u or v grids.
                ii = int(LOrivs_df.loc[rn,'col_py'])
                jj = int(LOrivs_df.loc[rn,'row_py'])
                
                uv = LOrivs_df.loc[rn,'uv']
                isign = LOrivs_df.loc[rn,'isign']
                idir = LOrivs_df.loc[rn,'idir']

                if uv == 'u' and isign == 1:
                    # River source on W side of rho cell
                    ax.plot(lon_u[jj,ii], lat_u[jj,ii],'>r')
                    if first_label_LO:
                       ax.plot(lon[jj,ii+1], lat[jj,ii+1],color='darkorange',
                               marker='*', linestyle = 'None', label='pre-existing LiveOcean River')
                       first_label_LO = False
                    else:
                        ax.plot(lon[jj,ii+1], lat[jj,ii+1],color='darkorange', marker='*')
                if uv == 'u' and isign == -1:
                    # River source on E side of rho cell
                    ax.plot(lon_u[jj,ii], lat_u[jj,ii],'<r')
                    ax.plot(lon[jj,ii], lat[jj,ii],color='darkorange', marker='*')
                if uv == 'v' and isign == 1:
                    # River source on S side of rho cell
                    ax.plot(lon_v[jj,ii], lat_v[jj,ii],'^b')
                    ax.plot(lon[jj+1,ii], lat[jj+1,ii],color='darkorange', marker='*')
                if uv == 'v' and isign == -1:
                    # River source on N side of rho cell
                    ax.plot(lon_v[jj,ii], lat_v[jj,ii],'vb')
                    ax.plot(lon[jj,ii], lat[jj,ii],color='darkorange', marker='*') 

#################################################################################
#                              Plot point sources                               #
#################################################################################
        
        # draw tracks for point sources
        if source_type == 'wwtp':

            # get new source locations
            ps_lon = []
            ps_lat = []
            for i,ind in enumerate(rowcol_df['col_py']):
                # remove nans (if point source is outside of the grid)
                if math.isnan(ind):
                    ps_lon = ps_lon + [np.nan]
                    ps_lat = ps_lat + [np.nan]
                else:
                    ps_lon = ps_lon + [X[int(ind)]]
                    indy = rowcol_df['row_py'][i]
                    ps_lat = ps_lat + [Y[int(indy)]]

            # plot original sources
            ax.plot(source_ds.lon.values,source_ds.lat.values, linestyle='none', markeredgecolor ='k',
                marker='D', color='yellowgreen', label='Original lat/lon coordinates')
            
            # add placed location
            ax.plot(ps_lon,ps_lat, linestyle='none',
                    marker='o', color='purple', label='Placed location in LO')
            
            # draw track from original location to placed location
            for ii in range(len(ps_lon)-1):
                ax.plot([source_ds.lon.values[ii], ps_lon[ii]],[source_ds.lat.values[ii], ps_lat[ii]],
                        color='purple', linewidth=1)

            # add source labels
            snames = source_ds.name.values
            for ii,source in enumerate(snames):
                text_lon = ps_lon[ii]
                text_lat = ps_lat[ii]+0.003
                ax.text(text_lon, text_lat, source, color = 'purple',
                        fontsize=12, horizontalalignment='center')

        # finalize figure
        ax.legend(loc='lower left', fontsize = 12)
        plt.show()

#################################################################################
#                         Save source info in csv file                          #
#################################################################################

    out_rfn = Ldir['grid'] / output_fn
    print('\nCreating ' + str(out_rfn))
    # drop any sources that were outside of the domain, and had values padded with nans
    rowcol_df = rowcol_df.dropna()
    rowcol_df.to_csv(out_rfn)

    return

#################################################################################
#                          Run function to place sources                        #
#################################################################################

# Run placement algorithm
traps_placement('riv',Ldir)
traps_placement('wwtp',Ldir)