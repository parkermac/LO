"""
Utility code to convert excel files to pandas DataFrames, in order to speed up
processing (?).

"""

import pandas as pd
from lo_tools import Lfun
Ldir = Lfun.Lstart()
import sys

source = 'ecology'
in_dir0 = Ldir['data'] / 'obs' / source

for fn in list(in_dir0.glob('*.xlsx')):

    if fn.name == 'ParkerMacCreadyCoreStationInfoFeb2018.xlsx':
        df = pd.read_excel(fn, index_col='Station')
        df.to_pickle(in_dir0 / 'sta_df.p')

    if fn.name == 'Parker_2006-2017_Nutrients.xlsx':
        df = pd.read_excel(fn, sheet_name='2006-2017', parse_dates = ['Date'])
        df.to_pickle(in_dir0 / 'bottle_2006_2017.p')

    if fn.name == 'ParkerMacCready1999-2016CTDDataMay2018.xlsx':
        df = pd.read_excel(fn, sheet_name='1999-2016Finalized_CTDResults', parse_dates = ['Date'])
        df.to_pickle(in_dir0 / 'ctd_1999_2016.p')

    if fn.name == 'ParkerMacCready2017CTDDataFeb2018.xlsx':
        df = pd.read_excel(fn, sheet_name='2017Provisional_CTDResults', parse_dates = ['Date'])
        df.to_pickle(in_dir0 / 'ctd_2017.p')

    if fn.name == 'ParkerMacCready2018CTDDOMar2020.xlsx':
        df = pd.read_excel(fn, sheet_name='2018_CTDDOResults', parse_dates = ['Date'])
        df.to_pickle(in_dir0 / 'ctd_2018.p')

    if fn.name == 'ParkerMacCready2019CTDDataFeb2020.xlsx':
        df = pd.read_excel(fn, sheet_name='2019Provisional_CTDResults', parse_dates = ['Date'])
        df.to_pickle(in_dir0 / 'ctd_2019.p')
        
    print('\n'+fn.name)
    print(df.columns)


"""
Screen output (helpful for names):
    
ParkerMacCreadyCoreStationInfoFeb2018.xlsx
Index(['Station', 'Lat_NAD83 (deg / dec_min)', 'Long_NAD83 (deg / dec_min)',
       'Desig', 'Descrip', 'Basin', 'Max_Depth'],
      dtype='object')

Parker_2006-2017_Nutrients.xlsx
Index(['ResultID', 'Project', 'Date', 'LocalTime', 'UTCDateTime', 'Year',
       'Month', 'Station', 'Niskin', 'Niskin Cast Rep', 'NomDepth',
       'Depth_Matching', 'Sampling Depth', 'PO4(uM)D', 'QC PO4_Lab',
       'SiOH4(uM)D', 'QC SiOH4_Lab', 'NO3(uM)D', 'QC NO3_Lab', 'NO2(uM)D',
       'QC NO2_Lab', 'NH4(uM)D', 'QC NH4_Lab', 'CTD Cast Rep', 'Type',
       'Replicate?', 'LabReplicateNumber', 'SampleFieldReplicateNumber', 'QA',
       'Filename', 'Bottle Number', 'Comments', 'Analysis Date'],
      dtype='object')

ParkerMacCready1999-2016CTDDataMay2018.xlsx
Index(['ResultID', 'Station', 'Date', 'Depth', 'Year', 'Temp', 'Temp QC Code',
       'Salinity', 'Salinity QC Code', 'Density', 'Density QC Code',
       'Chla_adjusted', 'Chla_raw QC Code', 'DO_adjusted', 'DO_sat_adjusted',
       'DO QC Code', 'Xmiss_25cm', 'Xmiss QC Code', 'Turbidity',
       'Turbidity QC Code', 'Rep'],
      dtype='object')

ParkerMacCready2017CTDDataFeb2018.xlsx
Index(['ResultID', 'Station', 'Date', 'Depth', 'Year', 'Temp', 'Salinity',
       'Density', 'Chla_adjusted', 'DO_raw', 'DO_sat', 'Xmiss_25cm',
       'Turbidity', 'Rep'],
      dtype='object')

ParkerMacCready2018CTDDOMar2020.xlsx
Index(['Station', 'Date', 'Depth', 'Rep', 'Temp', 'TempQC', 'TempQF', 'TempQA',
       'Salinity', 'SalinityQC', 'SalinityQF', 'SalinityQA', 'Density',
       'DensityQC', 'DensityQF', 'DensityQA', 'DO_raw', 'DO_adjusted',
       'DO_rawQC', 'DO_rawQF', 'DO_rawQA', 'DO_sat', 'DO_sat_adjusted',
       'DO_satQC', 'DO_satQF', 'DO_satQA'],
      dtype='object')

ParkerMacCready2019CTDDataFeb2020.xlsx
Index(['Station', 'Date', 'Depth', 'Rep', 'Temp', 'TempQC', 'TempQF', 'TempQA',
       'Salinity', 'SalinityQC', 'SalinityQF', 'SalinityQA', 'Density',
       'DensityQC', 'DensityQF', 'DensityQA', 'DO_raw', 'DO_rawQC', 'DO_rawQF',
       'DO_rawQA', 'DO_sat', 'DO_satQC', 'DO_satQF', 'DO_satQA', 'Fluor_adj',
       'Fluor_adjQC', 'Fluor_adjQF', 'Fluor_adjQA', 'FluorAdjComment',
       'Turbidity', 'TurbidityQC', 'TurbidityQF', 'TurbidityQA', 'Xmiss_25cm',
       'Xmiss_25cmQC', 'Xmiss_25cmQF', 'Xmiss_25cmQA'],
      dtype='object')
    NOTE: this spreadsheet also has a sheet with 2019 nutrient data!
"""