"""
Functions to pull river data from websites.

The input 'rs' is a pandas Series pulled from the river_info DataFrame.
The input 'days' is a tuple of datetimes.
"""

import sys
import pandas as pd
import urllib
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import requests
import bs4
import numpy as np

timeout = 10
print_exception = False
    
def get_usgs_data_custom(rs, days, temperature=False):
    """
    A wrapper for get_usgs_data_sub() that deals with any rivers
    requiring special handling.
    """
    if rs.name == 'skokomish':
        print(' combining to form Skokomish River '.center(60,'+'))
        rs.usgs = 12060500
        rs.ratio = 1.4417
        rs, qt1 = get_usgs_data(rs, days, temperature=temperature)
        got1 = rs.got_data
        rs.usgs = 12059500
        rs.ratio = 1.0
        rs, qt2 = get_usgs_data(rs, days, temperature=temperature)
        rs.got_data = rs.got_data and got1
        if rs.got_data:
            qt = qt1 + qt2
        else:
            qt = ''
    elif rs.name == 'hamma':
        print(' combining to form Hamma Hamma River '.center(60,'+'))
        rs.usgs = 12060500
        rs.ratio = 1.4417
        rs, qt1 = get_usgs_data(rs, days, temperature=temperature)
        got1 = rs.got_data
        rs.usgs = 12059500
        rs.ratio = 1.0
        rs, qt2 = get_usgs_data(rs, days, temperature=temperature)
        rs.got_data = rs.got_data and got1
        if rs.got_data:
            qt = 0.4125 * (qt1 + qt2)
        else:
            qt = ''
    else:
        print(' River not supported for custom extraction! '.center(60,'*'))
        qt = ''
        rs['got_data'] = False
    return rs, qt

def get_usgs_data(rs, days, temperature=False):
    """
    This gets USGS data for a past time period specified by 'days',
    a tuple of datetimes.  If 'days is empty then we get the most recent 6 days.
    The default is that it returns a pandas Series that is a timeseries of flow
    (m3/s).
    
    Use temperature=True to get temperature (degC) instead of flow.
    """
    # Set the time period to get.
    if len(days) == 2:
        time_str = ('&startDT=' + days[0].strftime('%Y-%m-%d')
            +'&endDT=' + days[1].strftime('%Y-%m-%d'))
    else:
        # This gets the most recent 6 days (daily intervals?)
        time_str = '&period=P6D'
    # Form the url.
    if temperature:
        url_str = ('http://waterservices.usgs.gov/nwis/dv/'
            + '?format=waterml,1.1&sites=' + str(int(rs.usgs))
            + time_str + '&parameterCd=00010')
    else:
        url_str = ('http://waterservices.usgs.gov/nwis/dv/'
            + '?format=waterml,1.1&sites=' + str(int(rs.usgs))
            + time_str + '&parameterCd=00060')
    try:
        # Get the XML.
        file = urllib.request.urlopen(url_str, timeout=timeout)
        tree = ET.parse(file)
        root = tree.getroot()
        # Extract the data from the XML.
        flag = True
        # aa = '{http://www.cuahsi.org/waterML/1.1/}'
        rt = root.tag
        aa = rt[rt.find('{'): rt.find('}') + 1]
        Q = []
        T = []
        for e0 in root.findall(".//"):
            if e0.tag == aa+'value':
                Q.append(float(e0.text))
                T.append(pd.to_datetime(e0.get('dateTime')))
            if e0.tag == aa+'unitCode' and flag:
                if temperature:
                    rs.temp_units = e0.text
                else:
                    rs.flow_units = e0.text
                flag = False
        qt = pd.Series(Q, index=T)
        rs, qt = fix_units(rs, qt)
        if not temperature:
            qt = float(rs.ratio) * qt
        # Note: this resampling of daily data just moves the timestamp to noon
        # of the day it started at.  Data is unchanged.
        qt = qt.resample('D', label='right', offset='-12h').mean()
        rs['got_data'] = True
    except Exception as e:
        qt = ''
        if print_exception:
            print(e)
        rs['got_data'] = False
    return rs, qt

def get_nws_data(rs):
    """
    This gets NWS forecast data.
    """
    url_str = ('http://www.nwrfc.noaa.gov/xml/xml.cgi?id=' + rs.nws +
               '&pe=HG&dtype=b&numdays=10')
    try:
        # get the XML
        file = urllib.request.urlopen(url_str, timeout=timeout)
        tree = ET.parse(file)
        root = tree.getroot()
        flag = True
        # NOTE: you find this tag by looking at any instance of e0.tag
        #aa = '{http://www.nwrfc.noaa.gov/xml/schemas/2004/03/hydromet_data}'
        rt = root.tag
        aa = rt[rt.find('{'): rt.find('}') + 1]
        Q = []
        T = []
        for e0 in root.findall(".//"):
            if e0.tag == aa+'observedData' or e0.tag == aa+'forecastData':
                for e in e0:
                    if e.tag == aa+'observedValue' or e.tag == aa+'forecastValue':
                        for ee in e:
                            if ee.tag == aa+'discharge':
                                Q.append(float(ee.text))
                                if flag:
                                    rs.flow_units = ee.get('units')
                                    # print(rs.flow_units)
                                    flag = False
                            if ee.tag == aa+'dataDateTime':
                                # print(ee.text)
                                T.append(pd.to_datetime(ee.text).tz_convert(None))
        qt = pd.Series(Q, index=T)
        rs, qt = fix_units(rs, qt)
        qt = float(rs.ratio) * qt
        qt = qt.resample('D', label='right', offset='-12h').mean()
        rs['got_data'] = True
    except Exception as e:
        qt = ''
        if print_exception:
            print(e)
        rs['got_data'] = False
    return rs, qt

def get_ec_data(rs, days, temperature=False):
    """
    Gets Environment Canada data, using code cribbed from:
    https://bitbucket.org/douglatornell/ecget/src/.
    NOTE: this will get data up through today, but can only go back 18 months into the past.
    To get longer records use get_ec_data_historical().
    """
    try:
        PARAM_IDS = {'discharge': 47, 'temperature': 5}
        if temperature:
            prm1 = PARAM_IDS['temperature']
        else:
            prm1 = PARAM_IDS['discharge']
        params = {
            'mode': 'Table',
            'type': 'realTime',
            'prm1': prm1,
            'prm2': -1,
            'stn': rs.ec,
            'startDate': days[0].strftime('%Y-%m-%d'),
            'endDate': days[1].strftime('%Y-%m-%d'),
        }
        DATA_URL = 'http://wateroffice.ec.gc.ca/report/real_time_e.html'
        DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
        response = requests.get(DATA_URL, params=params,
                                cookies=DISCLAIMER_COOKIE, timeout=timeout)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        table = soup.find('table')
        table_body = table.find('tbody')
        rows = table_body.find_all('tr')
        data = []
        for row in rows:
            cols = row.find_all('td')
            cols = [ele.text.strip() for ele in cols]
            data.append([ele for ele in cols if ele]) # Get rid of empty values
        d_dict = dict()
        for item in data:
            d_dict[pd.to_datetime(item[0])] = float(item[1].replace(',',''))
        qt = pd.Series(d_dict)
        rs, qt = fix_units(rs, qt)
        if not temperature:
            qt = float(rs.ratio) * qt
        qt = qt.resample('D', label='right', offset='-12h').mean()
        
        # NEW 2019.03.20 to deal with the problem that when you request date from
        # before 18 months ago it gives the most recent data instead.
        dt0_actual = qt.index[0]
        dt0_requested = days[0]
        import numpy as np
        if np.abs((dt0_actual - dt0_requested).days) >= 1:
            print('That date range was not available')
            rs['got_data'] = False
            qt = ''
        else:
            rs['got_data'] = True
    except Exception as e:
        qt = ''
        if print_exception:
            print(e)
        rs['got_data'] = False
    return rs, qt

def get_ec_data_historical(rs, year):
    # gets historical Environment Canada data, a year at a time,
    # using code cribbed from:
    # https://bitbucket.org/douglatornell/ecget/src/
    try:
        params = {
            'mode': 'Table',
            'type': 'h2oArc',
            'stn': rs.ec,
            'dataType': 'Daily',
            'parameterType': 'Flow',
            'year': str(year),
            'y1Max': '1',
            'y1Min': '1',
        }
        DATA_URL = 'http://wateroffice.ec.gc.ca/report/historical_e.html'
        DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
        response = requests.get(DATA_URL, params=params,
                                cookies=DISCLAIMER_COOKIE, timeout=timeout)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        if (str(year) + ' Daily Discharge') in soup.text:
            # we do the test above because the website will return the
            # table for the most recent year if the requested year is
            # missing
            table = soup.find('table')
            table_body = table.find('tbody')
            rows = table_body.find_all('tr')
            d_dict = dict()
            # Note that the sequence of the read loops below means that the
            # data are not in chronological order, so you need to sort later.
            for row in rows:
                # what day is it
                for ele in row.find_all('th'):
                    day = int(ele.text.strip())
                cols = row.find_all('td')
                cols = [ele.text.strip() for ele in cols] # a list of strings
                imo = 1
                for item in cols:
                    if (len(item) == 0) or (item=='-'):
                       pass
                    else:
                       this_data = float(item.split()[0].replace(',',''))
                       # the split call is to remove trailing letters
                       # that occasionally appear after the data
                       d_dict[datetime(year, imo, day, 12, 0, 0)] = this_data
                    imo += 1
            qt = pd.Series(d_dict)
            #rs.fix_units() # not needed
            qt = float(rs.ratio) * qt
            rs['got_data'] = True
        else:
            qt = ''
            rs['got_data'] = False
    except Exception as e:
        qt = ''
        if print_exception:
            print(e)
        rs['got_data'] = False
    return rs, qt

def fix_units(rs, qt):
    # fix flow units
    if rs.flow_units == 'kcfs':
        qt = qt*28.3168466
        rs.flow_units = 'm3/s'
    elif (rs.flow_units == 'cubic feet per second') or (rs.flow_units == 'ft3/s'):
        qt = qt*0.0283168466
        rs.flow_units = 'm3/s'
    elif rs.flow_units == 'm3/s':
        pass
    else:
        print('ERROR: Unrecognized flow units!')
        sys.exit()
    # fix temperature units
    if rs.temp_units == 'deg C':
        rs.temp_units = 'degC'
    elif rs.temp_units == 'degC':
        pass
    else:
        print('ERROR: Unrecognized temperature units!')
        sys.exit()
    return rs, qt
    
# testing
if __name__ == '__main__':
    print_exception = True
    import Lfun
    Ldir = Lfun.Lstart()
    # Load a dataframe with info for rivers to get
    ri_fn = Ldir['LOo'] / 'pre' / 'river' / 'cas6_v3' / 'river_info.csv'
    df = pd.read_csv(ri_fn, index_col='rname')
    
    if False:
        # test usgs standard case
        rn = 'skagit'
        print('\n'+(' testing usgs ' + rn).center(60,'-'))
        rs = df.loc[rn].copy()
        dsf = Ldir['ds_fmt']
        days = (datetime(2018,1,1), datetime(2018,1,10))
        rs, qt = get_usgs_data(rs, days)
        print(rs)
        print(qt)
        
    if False:
        # test usgs standard case for a forecast time period
        rn = 'columbia'
        print('\n'+(' testing usgs ' + rn).center(60,'-'))
        rs = df.loc[rn].copy()
        ds0 = datetime.now().strftime(Lfun.ds_fmt)
        dt0 = datetime.strptime(ds0,Lfun.ds_fmt) - timedelta(days=2.5)
        dt1 = datetime.strptime(ds0,Lfun.ds_fmt) + timedelta(days=4.5)
        days = (dt0, dt1)
        rs, qt = get_usgs_data(rs, days)
        print(rs)
        print(qt)
        
    if False:
        # test usgs temperature
        rn = 'skagit'
        print('\n'+(' testing usgs temperature ' + rn).center(60,'-'))
        rs = df.loc[rn].copy()
        days = (datetime(2018,1,1), datetime(2018,1,10))
        rs, qt = get_usgs_data(rs, days, temperature=True)
        print(rs)
        print(qt)
        
    if False:
        # test usgs for a "custom" river
        rn = 'hamma'
        rs = df.loc[rn].copy()
        days = (datetime(2021,4,16), datetime(2021,4,23))
        print('\n'+(' testing usgs custom ' + rn).center(60,'-'))
        rs = df.loc[rn].copy()
        rs, qt = get_usgs_data_custom(rs, days)
        print(rs)
        print(qt)
    
    if False:
        # test usgs for a river that is not on the custom list
        # should throw an error
        rn = 'skagit'
        rs = df.loc[rn].copy()
        days = (datetime(2018,1,1), datetime(2018,1,10))
        print('\n'+(' testing usgs custom ' + rn).center(60,'-'))
        rs = df.loc[rn].copy()
        rs, qt = get_usgs_data_custom(rs, days)
        print(rs)
        print(qt)
    
    if False:
        # test nws
        rn = 'puyallup'
        print('\n'+(' testing nws ' + rn).center(60,'-'))
        rs = df.loc[rn].copy()
        rs, qt = get_nws_data(rs)
        print(rs)
        print(qt)
        
    if False:
        # test ec
        rn = 'fraser'
        for year in range(2021,2023):
            print('\n'+(str(year) + ': testing ec ' + rn).center(60,'-'))
            for mo in [1,3,5,7,9,11]:
                rs = df.loc[rn].copy()
                days = (datetime(year,mo,1), datetime(year,mo,10))
                rs, qt = get_ec_data(rs, days)
                print(rs)
                print(qt)
    
    if False:
        # test ec temperature
        rn = 'fraser'
        for year in range(2018,2022):
            print('\n'+(str(year) + ': testing ec temperature ' + rn).center(60,'-'))
            rs = df.loc[rn].copy()
            days = (datetime(year,1,1), datetime(year,1,10))
            rs, qt = get_ec_data(rs, days, temperature=True)
            print(rs)
            print(qt)
    
    if False:
        # test ec historical
        rn = 'fraser'
        for year in range(2018,2022):
            print('\n'+(str(year) + ': testing ec historical ' + rn).center(60,'-'))
            rs = df.loc[rn].copy()
            rs, qt = get_ec_data_historical(rs, year)
            print(rs)
            print(qt)
        
    
    
    
