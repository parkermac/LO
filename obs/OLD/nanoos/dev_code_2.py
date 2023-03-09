"""
Code to help with the development of the nanoos processing.

"""

import pandas as pd
import numpy as np
import gsw

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import bot_fun
v_dict = bot_fun.v_dict

# BOTTLE
source = 'nanoos'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source
year_list = [2017, 2021]

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

a = []
for year in year_list:
    ys = str(year)
    
    # name output files
    out_fn = out_dir / (str(year) + '.p')
    info_out_fn = out_dir / ('info_' + str(year) + '.p')
    
    in_dirs = list(in_dir0.glob('*'+ys+'*'))
    fn_list = []
    for in_dir in in_dirs:
        fn_list.append(list(in_dir.glob('*labupcast*'))[0])
        
    cid0 = 0
    DF = pd.DataFrame()
    for fn in fn_list:
        print('\n'+fn.name)
        # Initial reading of the files.
        if fn.name == 'SalishCruise_April_2021_labupcast.xlsx':
            # Hack required to get this cruise to read correctly, even after hand editing
            # to remove a lot of bad lines at the end of the excel file.
            df0 = pd.read_excel(fn, dtype={'DATE_UTC':str, 'TIME_UTC':str})
            # make a time column by hand
            dstr = df0.DATE_UTC
            tstr = df0.TIME_UTC
            tlist = []
            for ii in range(len(dstr)):
                try:
                    tlist.append(pd.to_datetime(dstr[ii][:10] + ' ' + tstr[ii][:8]))
                except:
                    tlist.append(np.nan)
            df0['time'] = tlist
        else:
            try:
                df0 = pd.read_excel(fn, parse_dates={'time':['DATE_UTC', 'TIME_UTC']})
            except ValueError:
                df0 = pd.read_excel(fn, parse_dates={'time':['Date_UTC', 'Time_UTC']})
        df0 = df0.dropna(axis=0, how='all') # drop rows with no good data
        df0 = df0[df0.time.notna()] # drop rows with bad time
        df0 = df0.reset_index()
        # Create a DataFrame with only known variables (those which have non-empty values
        # in bot_fun.v_dict).
        df = pd.DataFrame()
        for v in df0.columns:
            if v in v_dict.keys():
                if len(v_dict[v]) > 0:
                    df[v_dict[v]] = df0[v]
        df['time'] = df0['time']
        
        # Now proceed with the processing to get a single DataFrame for the year.
        
        # add the "cid" (cast ID) column
        #
        # Note that we will save the field "name" for station number, since this dataset has
        # repeat stations which is helpful for plotting sections. Then we will generate our own
        # cid, a unique one for each cast, being careful to keep them unique for the collection
        # of cruises in this year, even though a station may be repeated on all cruises.
        #
        # We will also save the field "cruise" as a convenient way to select a collection of
        # casts from a given month.
        df['cid'] = np.nan
        cid = cid0
        for name in df.name.unique():
            df.loc[df.name==name,'cid'] = cid
            cid += 1
        for cid in df.cid.unique():
            # Check that there are not two different casts associated with the same station
            # by looking for large time differences. Pretty ad hoc.
            time_diff = df[df.cid==cid].time.values[-1] - df[df.cid==cid].time.values[0]
            time_diff = pd.to_timedelta(time_diff)
            if time_diff.days > 1 or time_diff.days < -1:
                print('Station %d has time diff of %d days' % (cid, time_diff.days))
            # RESULT: the time_diffs are all zero, so it appears that in this database
            # the Station field is a unique cast identifier.

            # Force certain fields to be the same throughout the cast.
            df.loc[df.cid==cid,'lon'] = - df[df.cid==cid].lon.values[0] # also fix sign
            df.loc[df.cid==cid,'lat'] = df[df.cid==cid].lat.values[0]
            df.loc[df.cid==cid,'time'] = df[df.cid==cid].time.values[0]
                
        # Next make derived quantities and do unit conversions
    
        # (1) Create CT, SA, and z
        # - pull out variables
        SP = df.SP.to_numpy()
        pt = df.PT.to_numpy()
        p = df['P (dbar)'].to_numpy()
        lon = -(np.abs(df.lon.to_numpy())) # force sign to be negative
        df['lon'] = lon
        lat = df.lat.to_numpy()
        # - do the conversions
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        CT = gsw.CT_from_pt(SA, pt)
        z = gsw.z_from_p(p, lat)
        # - add the results to the DataFrame
        df['SA'] = SA
        df['CT'] = CT
        df['z'] = z
    
        # (2) units
        if 'DO (mg/L)' in df.columns:
            df['DO (uM)'] = (1000/32) * df['DO (mg/L)']
        elif 'DO (umol/kg)' in df.columns:
            rho = gsw.rho(SA,CT,p)
            df['DO (uM)'] = (rho/1000) * df['DO (umol/kg)']
            
        # Screen output to check on the results on the results.
        print(df.columns)
        print(df.loc[0,'time'])
        print(df.loc[len(df)-1,'time'])
        
    
        # (3) retain only selected variables
        cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
            'CT', 'SA', 'DO (uM)',
            'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
            'ChlA (ug/L)', 'TA (umol/kg)', 'DIC (umol/kg)']
        this_cols = [item for item in cols if item in df.columns]
        df = df[this_cols]
        DF = pd.concat([DF,df])
            
        print(' - processed %d casts' % ( len(df.cid.unique()) ))
        cid0 = df.cid.max() + 1
            
        # Sort the result by time, and sort each cast to be bottom to top
        DF = DF.sort_values(['time','z'], ignore_index=True)
        
        # Rework cid to also be increasing in time
        a = DF[['time','cid']].copy()
        a['cid_alt'] = np.nan
        ii = 0
        for t in a.time.unique():
            a.loc[a.time==t,'cid_alt'] = ii
            ii += 1
        DF['cid'] = a['cid_alt'].copy()

        # clean up spaces in cruise names
        a = DF['cruise'].to_list()
        aa = [item.strip() for item in a]
        DF['cruise'] = aa
        
        # Hack: fix a location typo
        if year == 2021:
            DF.loc[(DF.cruise=='RC0051') & (DF.name==5), 'lat'] = 47.883

        # Save the data
        DF.to_pickle(out_fn)
    
        # Also pull out a dateframe with station info to use for model cast extractions.
        ind = DF.cid.unique()
        col_list = ['lon','lat','time','name','cruise']
        info_df = pd.DataFrame(index=ind, columns=col_list)
        for cid in DF.cid.unique():
            info_df.loc[cid,col_list] = DF.loc[DF.cid==cid,col_list].iloc[0,:]
        info_df.index.name = 'cid'
        info_df['time'] = pd.to_datetime(info_df['time'])
        info_df.to_pickle(info_out_fn)
