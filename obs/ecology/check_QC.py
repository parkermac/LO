"""
Utility code to check on QC codes in the Ecology data.

The code is 2 for good data, but it is in various formats which
makes processing complicated. In some cases it is a string like
'2_0_2' and in others it is the integer 2, and in others there
is a nan, probably because it was missing in the original excel
file.

We work through all dataframes, putting in nan's where the flag is
anything but 2, and saving it as a new pickled dataframe []_fixed.p.
"""

import pandas as pd
import numpy as np
from lo_tools import Lfun
Ldir = Lfun.Lstart()
import sys

source = 'ecology'
in_dir0 = Ldir['data'] / 'obs' / source

fn_list = list(in_dir0.glob('*.p'))
fn_list = [item for item in fn_list if 'fixed' not in item.name]
fn_list = [item for item in fn_list if 'sta_df' not in item.name]

for fn in fn_list:

    df = pd.read_pickle(fn)
        
    c = [item for item in df.columns if 'QC' in item]
    
    if len(c) > 0:
        print('\n'+fn.name)
        print(df.columns)
        for cc in c:
            a = df[cc].to_list()
            x = []
            for aa in a:
                if isinstance(aa,(int,float)):
                    if aa == 2:
                        x.append(1)
                    else:
                        x.append(0)
                elif isinstance(aa,str):
                    if aa[0] == '2':
                        x.append(1)
                    else:
                        x.append(0)
                        
            xx = np.array(x)
            ngood = xx.sum()
            nn = len(xx)
            nbad = nn - ngood
            # also attempt to automate finding the name of the
            # associated data column
            C = cc.replace('QC','').strip().split()[0].replace('_Lab','(uM)D')
            if C == 'DO':
                C = 'DO_adjusted' # hack for one exception
            if C == 'Chla_raw':
                C = 'Chla_adjusted' # hack for one exception
            if C == 'Xmiss':
                C = 'Xmiss_25cm' # hack for one exception
            if C not in df.columns:
                r = 'XX'
            else:
                r = 'OK'
            print('%2s %20s : %10s || %d bad out of %d' % (r, cc, C, nbad, nn))
            
            # implement masking and save
            df.loc[xx==0,C] = np.nan
            
    new_name = fn.name.replace('.p', '_fixed.p')
    df.to_pickle(in_dir0 / new_name)
    
"""
Screen output:

ctd_2019.p
Index(['Station', 'Date', 'Depth', 'Rep', 'Temp', 'TempQC', 'TempQF', 'TempQA',
       'Salinity', 'SalinityQC', 'SalinityQF', 'SalinityQA', 'Density',
       'DensityQC', 'DensityQF', 'DensityQA', 'DO_raw', 'DO_rawQC', 'DO_rawQF',
       'DO_rawQA', 'DO_sat', 'DO_satQC', 'DO_satQF', 'DO_satQA', 'Fluor_adj',
       'Fluor_adjQC', 'Fluor_adjQF', 'Fluor_adjQA', 'FluorAdjComment',
       'Turbidity', 'TurbidityQC', 'TurbidityQF', 'TurbidityQA', 'Xmiss_25cm',
       'Xmiss_25cmQC', 'Xmiss_25cmQF', 'Xmiss_25cmQA'],
      dtype='object')
OK               TempQC :       Temp || 78 bad out of 56440
OK           SalinityQC :   Salinity || 109 bad out of 56440
OK            DensityQC :    Density || 109 bad out of 56440
OK             DO_rawQC :     DO_raw || 541 bad out of 56440
OK             DO_satQC :     DO_sat || 541 bad out of 56440
OK          Fluor_adjQC :  Fluor_adj || 214 bad out of 56440
OK          TurbidityQC :  Turbidity || 215 bad out of 56440
OK         Xmiss_25cmQC : Xmiss_25cm || 44 bad out of 56440

ctd_1999_2016_fixed.p
Index(['ResultID', 'Station', 'Date', 'Depth', 'Year', 'Temp', 'Temp QC Code',
       'Salinity', 'Salinity QC Code', 'Density', 'Density QC Code',
       'Chla_adjusted', 'Chla_raw QC Code', 'DO_adjusted', 'DO_sat_adjusted',
       'DO QC Code', 'Xmiss_25cm', 'Xmiss QC Code', 'Turbidity',
       'Turbidity QC Code', 'Rep'],
      dtype='object')
OK         Temp QC Code :       Temp || 0 bad out of 855085
OK     Salinity QC Code :   Salinity || 6172 bad out of 855085
OK      Density QC Code :    Density || 6133 bad out of 855085
OK     Chla_raw QC Code : Chla_adjusted || 147326 bad out of 855085
OK           DO QC Code : DO_adjusted || 15478 bad out of 855085
OK        Xmiss QC Code : Xmiss_25cm || 19888 bad out of 855085
OK    Turbidity QC Code :  Turbidity || 349796 bad out of 855085

ctd_2018_fixed.p
Index(['Station', 'Date', 'Depth', 'Rep', 'Temp', 'TempQC', 'TempQF', 'TempQA',
       'Salinity', 'SalinityQC', 'SalinityQF', 'SalinityQA', 'Density',
       'DensityQC', 'DensityQF', 'DensityQA', 'DO_raw', 'DO_adjusted',
       'DO_rawQC', 'DO_rawQF', 'DO_rawQA', 'DO_sat', 'DO_sat_adjusted',
       'DO_satQC', 'DO_satQF', 'DO_satQA'],
      dtype='object')
OK               TempQC :       Temp || 213 bad out of 59041
OK           SalinityQC :   Salinity || 312 bad out of 59041
OK            DensityQC :    Density || 312 bad out of 59041
OK             DO_rawQC :     DO_raw || 5655 bad out of 59041
OK             DO_satQC :     DO_sat || 5655 bad out of 59041

ctd_1999_2016.p
Index(['ResultID', 'Station', 'Date', 'Depth', 'Year', 'Temp', 'Temp QC Code',
       'Salinity', 'Salinity QC Code', 'Density', 'Density QC Code',
       'Chla_adjusted', 'Chla_raw QC Code', 'DO_adjusted', 'DO_sat_adjusted',
       'DO QC Code', 'Xmiss_25cm', 'Xmiss QC Code', 'Turbidity',
       'Turbidity QC Code', 'Rep'],
      dtype='object')
OK         Temp QC Code :       Temp || 0 bad out of 855085
OK     Salinity QC Code :   Salinity || 6172 bad out of 855085
OK      Density QC Code :    Density || 6133 bad out of 855085
OK     Chla_raw QC Code : Chla_adjusted || 147326 bad out of 855085
OK           DO QC Code : DO_adjusted || 15478 bad out of 855085
OK        Xmiss QC Code : Xmiss_25cm || 19888 bad out of 855085
OK    Turbidity QC Code :  Turbidity || 349796 bad out of 855085

ctd_2019_fixed.p
Index(['Station', 'Date', 'Depth', 'Rep', 'Temp', 'TempQC', 'TempQF', 'TempQA',
       'Salinity', 'SalinityQC', 'SalinityQF', 'SalinityQA', 'Density',
       'DensityQC', 'DensityQF', 'DensityQA', 'DO_raw', 'DO_rawQC', 'DO_rawQF',
       'DO_rawQA', 'DO_sat', 'DO_satQC', 'DO_satQF', 'DO_satQA', 'Fluor_adj',
       'Fluor_adjQC', 'Fluor_adjQF', 'Fluor_adjQA', 'FluorAdjComment',
       'Turbidity', 'TurbidityQC', 'TurbidityQF', 'TurbidityQA', 'Xmiss_25cm',
       'Xmiss_25cmQC', 'Xmiss_25cmQF', 'Xmiss_25cmQA'],
      dtype='object')
OK               TempQC :       Temp || 78 bad out of 56440
OK           SalinityQC :   Salinity || 109 bad out of 56440
OK            DensityQC :    Density || 109 bad out of 56440
OK             DO_rawQC :     DO_raw || 541 bad out of 56440
OK             DO_satQC :     DO_sat || 541 bad out of 56440
OK          Fluor_adjQC :  Fluor_adj || 214 bad out of 56440
OK          TurbidityQC :  Turbidity || 215 bad out of 56440
OK         Xmiss_25cmQC : Xmiss_25cm || 44 bad out of 56440

bottle_2006_2017.p
Index(['ResultID', 'Project', 'Date', 'LocalTime', 'UTCDateTime', 'Year',
       'Month', 'Station', 'Niskin', 'Niskin Cast Rep', 'NomDepth',
       'Depth_Matching', 'Sampling Depth', 'PO4(uM)D', 'QC PO4_Lab',
       'SiOH4(uM)D', 'QC SiOH4_Lab', 'NO3(uM)D', 'QC NO3_Lab', 'NO2(uM)D',
       'QC NO2_Lab', 'NH4(uM)D', 'QC NH4_Lab', 'CTD Cast Rep', 'Type',
       'Replicate?', 'LabReplicateNumber', 'SampleFieldReplicateNumber', 'QA',
       'Filename', 'Bottle Number', 'Comments', 'Analysis Date'],
      dtype='object')
OK           QC PO4_Lab :   PO4(uM)D || 0 bad out of 11809
OK         QC SiOH4_Lab : SiOH4(uM)D || 1 bad out of 11809
OK           QC NO3_Lab :   NO3(uM)D || 4 bad out of 11809
OK           QC NO2_Lab :   NO2(uM)D || 0 bad out of 11809
OK           QC NH4_Lab :   NH4(uM)D || 1 bad out of 11809

bottle_2006_2017_fixed.p
Index(['ResultID', 'Project', 'Date', 'LocalTime', 'UTCDateTime', 'Year',
       'Month', 'Station', 'Niskin', 'Niskin Cast Rep', 'NomDepth',
       'Depth_Matching', 'Sampling Depth', 'PO4(uM)D', 'QC PO4_Lab',
       'SiOH4(uM)D', 'QC SiOH4_Lab', 'NO3(uM)D', 'QC NO3_Lab', 'NO2(uM)D',
       'QC NO2_Lab', 'NH4(uM)D', 'QC NH4_Lab', 'CTD Cast Rep', 'Type',
       'Replicate?', 'LabReplicateNumber', 'SampleFieldReplicateNumber', 'QA',
       'Filename', 'Bottle Number', 'Comments', 'Analysis Date'],
      dtype='object')
OK           QC PO4_Lab :   PO4(uM)D || 0 bad out of 11809
OK         QC SiOH4_Lab : SiOH4(uM)D || 1 bad out of 11809
OK           QC NO3_Lab :   NO3(uM)D || 4 bad out of 11809
OK           QC NO2_Lab :   NO2(uM)D || 0 bad out of 11809
OK           QC NH4_Lab :   NH4(uM)D || 1 bad out of 11809

ctd_2018.p
Index(['Station', 'Date', 'Depth', 'Rep', 'Temp', 'TempQC', 'TempQF', 'TempQA',
       'Salinity', 'SalinityQC', 'SalinityQF', 'SalinityQA', 'Density',
       'DensityQC', 'DensityQF', 'DensityQA', 'DO_raw', 'DO_adjusted',
       'DO_rawQC', 'DO_rawQF', 'DO_rawQA', 'DO_sat', 'DO_sat_adjusted',
       'DO_satQC', 'DO_satQF', 'DO_satQA'],
      dtype='object')
OK               TempQC :       Temp || 213 bad out of 59041
OK           SalinityQC :   Salinity || 312 bad out of 59041
OK            DensityQC :    Density || 312 bad out of 59041
OK             DO_rawQC :     DO_raw || 5655 bad out of 59041
OK             DO_satQC :     DO_sat || 5655 bad out of 59041    
"""

        
        
