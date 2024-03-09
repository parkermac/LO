"""
Code to re-process LineP CTD data. It is not compiled, you have to download year by year first

NOTE: they are two batch files, we need to handle them separately

(1) between 1993 and 2006, the 'date' and 'time' was not written as data columns in each .ctd file (ASCII format)
same for 'lat','lon', 'cruise', and 'station' and 'cast',
so we need to extract these info from a given file and append them as data columns, and save them to .csv.

(2) after 2006, ctd files are .csv, which is easier to handle

(3) there are multiple files for a given year, we want to combine those files to yearly

"""

import pandas as pd
import glob
import os
from lo_tools import Lfun
Ldir = Lfun.Lstart()

# CTD
source = 'LineP'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source /'ctd_uncompiled'
# output location
out_dir0 = Ldir['data'] / 'obs' / source / otype
Lfun.make_dir(out_dir0)

testing = False

if testing:
    year_list = [1993]
    # year_list = list(range(1993,2006+1))
    # year_list = list(range(2007,2023+1))
else:
    year_list = list(range(1993,2023+1))

ls = os.listdir(in_dir0)
vn_list = ['Pressure', 'Temperature', 'Salinity','Transmissivity','Oxygen']

for year in year_list:
    ys = str(year)
    print('\n'+ys)
    if year <= 2006:
        df = pd.DataFrame()
        for folder in ls:
            if ys in folder:
                in_dir = os.path.join(in_dir0, folder)
                print(in_dir)
                file_list = os.listdir(in_dir)
                file_list= sorted(file_list, key=lambda x: int(x.split('.')[0].split('-')[2]))
                for file in file_list:
                    # print(file)
                    in_file = os.path.join(in_dir, file)
                    print(in_file)
                    """ 
                    these .ctd files have '*END OF HEADER' before the data variables, 
                    so I created skiprows for each file 
                    """
                    with open(in_file, 'r', encoding='latin1') as file:
                        header_info = {}
                        for line_num, line in enumerate(file):
                            if line.strip() == '*END OF HEADER':
                                skiprows = line_num + 1  # Skip the line with *END OF HEADER
                                break
                            if line.startswith('    START TIME'):
                                start_time = line.split('UTC')[1].strip()
                                DateTime = start_time.split(',')[0].strip().split()
                                header_info['Date [UTC]'] = DateTime[0].strip()
                                header_info['Time [UTC]'] = DateTime[1].strip()
                            elif line.startswith('    MISSION'):
                                header_info['Mission'] = line.split(':')[1].strip()
                            elif line.startswith('    PROJECT'):
                                header_info['Project'] = line.split(':')[1].strip()
                            elif line.startswith('    SCIENTIST'):
                                header_info['Scientist'] = line.split(':')[1].strip()
                            elif line.startswith('    PLATFORM'):
                                header_info['Platform'] = line.split(':')[1].strip()
                            elif line.startswith('    STATION'):
                                header_info['Station'] = line.split(':')[1].strip()
                            elif line.startswith('    EVENT NUMBER'):
                                header_info['Event_number'] = line.split(':')[1].strip()
                            elif line.startswith('    LATITUDE'):
                                # split degrees and minutes
                                lat_deg, lat_min = line.split(':')[1].split('N')[0].strip().split()
                                header_info['LATITUDE_deg'] = lat_deg
                                header_info['LATITUDE_min'] = lat_min
                            elif line.startswith('    LONGITUDE'):
                                # split degrees and minutes
                                lon_deg, lon_min = line.split(':')[1].split('W')[0].strip().split()
                                header_info['LONGITUDE_deg'] = lon_deg
                                header_info['LONGITUDE_min'] = lon_min

                    print(f'skiprows = {skiprows}')
                    df0 = pd.read_csv(in_file, delimiter=r'\s+|,', skiprows=skiprows, engine='python', names=vn_list, encoding='latin1')
                    # drop rows with all NaN data
                    df0 = df0.dropna(axis=0, how='all')
                    df0 = df0.reset_index(drop=True)

                    # Add header information as columns to the DataFrame
                    for key, value in header_info.items():
                        df0[key] = value

                    df = pd.concat([df, df0], ignore_index=True)

                    # Calculate decimal Latitude
                    df['LATITUDE_deg'] = pd.to_numeric(df['LATITUDE_deg'])
                    df['LATITUDE_min'] = pd.to_numeric(df['LATITUDE_min'])
                    df['Latitude'] = df['LATITUDE_deg'] + df['LATITUDE_min'] / 60
                    # Calculate decimal Longitude
                    df['LONGITUDE_deg'] = pd.to_numeric(df['LONGITUDE_deg'])
                    df['LONGITUDE_min'] = pd.to_numeric(df['LONGITUDE_min'])
                    df['Longitude'] =  (df['LONGITUDE_deg'] + df['LONGITUDE_min'] / 60) * (-1)
        v_dict = {

            'Mission':'Mission',
            'Project':'Project',
            'Platform': 'Platform',
            'Scientist':'Scientist',
            'DateTime [UTC]':'',
            'Date [UTC]': 'Date [UTC]',
            'Time [UTC]': 'Time [UTC]',
            'Event_number':'Event_number',
            'Station':'Station',
            'Longitude':'Longitude',
            'Latitude':'Latitude',
            'Pressure':'Pressure [dbar]',
            'Temperature':'Temperature',
            'Salinity':'Salinity',
            'Transmissivity':'',
            'Oxygen':'Oxygen [umol/kg]'

        }

        df1 = pd.DataFrame()
        for v in v_dict.keys():
            if v in df.columns:
                if len(v_dict[v]) > 0:
                    df1[v_dict[v]] = df[v]

        df1.to_csv(out_dir0/f'{year}.csv', index=False)

    elif year > 2006:
        file_list = list(in_dir0.glob(f'*{year}*.csv'))
        file_list= sorted(file_list, key=lambda x: int(x.stem.split('-')[1]))
        print(file_list)

        """ 
        these .csv files have data starting from different rows, 
        so I manually looked through these files,
        and created this skiprows_list
        """
        skiprows_dict = {
            '2007':[11,13,13], #2007
            '2008':[12,15,15], #2008
            '2009':[42,22],    #2009
            '2010':[25,21,38], #2010
            '2011':[25,15,16], #2011
            '2012':[18,16,16], #2012
            '2013':[16,29,26], #2013
            '2014':[18,17,22], #2014
            '2015':[15,19,18], #2015
            '2016':[19,24,24], #2016
            '2017':[28,27,25], #2017
            '2018':[22,22,29], #2018
            '2019':[22,27,26], #2019
            '2020':[28,27,30], #2020
            '2021':[27,29,34], #2021
            '2022':[34,35],    #2022
            '2023':[48]       #2023
        }

        v_dict = {
            'File Name': '',
            'Zone': '',
            'FIL:START TIME YYYY/MM/DD':'Date [UTC]',
            ' HH:MM':'Time [UTC]',
            'LOC:EVENT_NUMBER':'Event_number',
            'LOC:LATITUDE':'Latitude',
            'LOC:LONGITUDE':'Longitude',
            'LOC:STATION':'Station',
            'INS:LOCATION':'',
            'Pressure:CTD [dbar]':'Pressure [dbar]',
            'Temperature:CTD [deg_C_(ITS90)]':'Temperature',
            'Salinity:CTD [PSS-78]':'Salinity',
            'Sigma-t:CTD [kg/m^3]':'Sigma-t [kg/m^3]',
            'Transmissivity:CTD:650 [*/m]':'',
            'Transmissivity:CTD:530 [*/m]':'',
            'Transmissivity:CTD [*/m]':'',
            'Oxygen:Dissolved:CTD:Volume [ml/l]':'Oxygen [ml/l]',
            'Oxygen:Dissolved:CTD:Mass [µmol/kg]':'Oxygen [umol/kg]',
            'Fluorescence:CTD:Seapoint [mg/m^3]':'Fluorescence [mg/m^3]',
            'PAR:CTD [µE/m^2/sec]':''
        }

        ii = 0
        df = pd.DataFrame()
        for file in file_list:
            print(file)
            if ys in skiprows_dict.keys():
                print(ii)
                df0 = pd.read_csv(file, skiprows=skiprows_dict[ys][ii], delimiter=',', engine='python', encoding='latin1')
                df1 = pd.DataFrame()
                for v in df0.columns:
                    if v in v_dict.keys():
                        if len(v_dict[v]) > 0:
                            df1[v_dict[v]] = df0[v]
                # drop rows with all NaN data
                df1 = df1.dropna(axis=0, how='all')
                df1 = df1.reset_index(drop=True)
                df = pd.concat([df, df1], ignore_index=True)

                ii +=1
                print('\n')

        df.to_csv(out_dir0/f'{year}.csv', index=False)

