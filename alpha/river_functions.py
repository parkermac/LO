"""
Functions to pull river data from websites.
"""

import pandas as pd
import urllib
import xml.etree.ElementTree as ET
from datetime import datetime

# "rs" is a pandas Series pulled from the river_info DataFrame
    
def get_usgs_data(rs, days):
    # this deals with rivers that need gauge combinations
    if rs.name == 'skokomish':
        print('+++ combining to form Skokomish River +++')
        rs.usgs = 12060500
        rs.ratio = 1.4417
        qt1 = get_usgs_data_sub(rs, days)
        rs.usgs = 12059500
        rs.ratio = 1.0
        qt2 = get_usgs_data_sub(rs, days)
        qt = qt1 + qt2
    elif rs.name == 'hamma':
        print('+++ combining to form Hamma Hamma River +++')
        rs.usgs = 12060500
        rs.ratio = 1.4417
        qt1 = get_usgs_data_sub(rs, days)
        rs.usgs = 12059500
        rs.ratio = 1.0
        qt2 = get_usgs_data_sub(rs, days)
        qt = 0.4125 * (qt1 + qt2)
    else:
        qt = get_usgs_data_sub(rs, days)
    return qt

def get_usgs_data_sub(rs, days):
    # This gets USGS data for a past time period specfied by
    # the tuple of datetimes "days".  If "days" is empty
    # then we get the most recent 6 days.

    # Set the time period to get.
    if len(days) == 2:
        time_str = ('&startDT=' + days[0].strftime('%Y-%m-%d')
            +'&endDT=' + days[1].strftime('%Y-%m-%d'))
    else:
        # This gets the most recent 6 days (daily intervals?)
        time_str = '&period=P6D'
    gage = int(rs.usgs)
    # Form the url.
    url_str = ('http://waterservices.usgs.gov/nwis/dv/'
        + '?format=waterml,1.1&sites=' + str(gage)
        + time_str + '&parameterCd=00060')
    try:
        # Get the XML.
        file = urllib.request.urlopen(url_str, timeout=10)
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
                rs.flow_units = e0.text
                flag = False
        qt = pd.Series(Q, index=T)
        rs, qt = fix_units(rs, qt)
        qt = float(rs.ratio) * qt
        # Note: this resampling of daily data just moves the timestamp to noon
        # of the day it started at.  Data is unchanged.
        qt = qt.resample('D', label='right', offset='-12h').mean()
        # rs.got_data = True
        print('success')
    except Exception as e:
        qt = ''
        print(e)
    
    return rs, qt

def get_nws_data(rs):
    # This gets NWS forecast data.
    url_str = ('http://www.nwrfc.noaa.gov/xml/xml.cgi?id=' +
               rs.nws +
               '&pe=HG&dtype=b&numdays=10')
    try:
        # get the XML
        file = urllib.request.urlopen(url_str, timeout=10)
        tree = ET.parse(file)
        root = tree.getroot()
        flag = True
        # NOTE: you find this tag by looking at any instance of e0.tag
        #aa = '{http://www.nwrfc.noaa.gov/xml/schemas/2004/03/hydromet_data}'
        rt = root.tag
        aa = rt[rt.find('{'): rt.find('}') + 1]
        for e0 in root.findall(".//"):
            if e0.tag == aa+'observedData' or e0.tag == aa+'forecastData':
                for e in e0:
                    if e.tag == aa+'observedValue' or e.tag == aa+'forecastValue':
                        for ee in e:
                            if ee.tag == aa+'discharge':
                                rs.Q.append(float(ee.text))
                                if flag:
                                    rs.flow_units = ee.get('units')
                                    flag = False
                            if ee.tag == aa+'dataDateTime':
                                rs.T.append(pd.to_datetime(ee.text))
        rs.qt = pd.Series(rs.Q, index=rs.T)
        rs.fix_units()
        rs.qt = float(rs.ratio) * rs.qt
        rs.qt = rs.qt.resample('D', label='right', loffset='-12h').mean()
        rs.got_data = True
        rs.memo = 'success'
    except:
        rs.memo = 'problem downloading or parsing data from XML'
    rs.memo = (rs.memo + ' NWS')

def get_ec_data(rs, days, timeout=10):
    # gets Environment Canada data, using code cribbed from:
    #https://bitbucket.org/douglatornell/ecget/src/
    # NOTE: this will get data up through today, but can only go back
    # 18 months into the past.
    # To get longer records use get_ec_data_historical.
    import requests
    import bs4
    try:
        PARAM_IDS = {'discharge': 47,}

        params = {
            'mode': 'Table',
            'type': 'realTime',
            'prm1': PARAM_IDS['discharge'],
            'prm2': -1,
            'stn': rs.ec,
            'startDate': days[0].strftime('%Y-%m-%d'),
            'endDate': days[1].strftime('%Y-%m-%d'),
        }
        DATA_URL = 'http://wateroffice.ec.gc.ca/report/real_time_e.html'
        DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
        response = requests.get(DATA_URL, params=params,
                                cookies=DISCLAIMER_COOKIE)
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
        rs.qt = pd.Series(d_dict)
        rs.flow_units = '$m^{3}s^{-1}$'
        #rs.fix_units() # not needed
        rs.qt = float(rs.ratio) * rs.qt
        rs.qt = rs.qt.resample('D', label='right', loffset='-12h').mean()
        # NEW 2019.03.20 to deal with the problem that when you request date from
        # before 18 months ago it gives the most recent data instead.
        dt0_actual = rs.qt.index[0]
        dt0_requested = days[0]
        import numpy as np
        if np.abs((dt0_actual - dt0_requested).days) >= 1:
            memo = 'That date range was not available'
            qt = ''
        else:
            rs.got_data = True
            rs.memo = 'success'
    except:
        rs.memo = 'problem parsing data from soup'
    rs.memo = (rs.memo + ' EC')

def get_ec_data_historical(rs, year):
    # gets Environment Canada data, using code cribbed from:
    #https://bitbucket.org/douglatornell/ecget/src/
    # NOTE: this will get data up through the end of 2016.
    import requests
    import bs4
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
                                cookies=DISCLAIMER_COOKIE)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        if (str(year) + ' Daily Discharge') in soup.text:
            # we do the test above because the website will return the
            # table for the most recent year if the requested year is
            # missing
            table = soup.find('table')
            table_body = table.find('tbody')
            rows = table_body.find_all('tr')
            d_dict = dict()
            for row in rows:
                # what day is it
                for ele in row.find_all('th'):
                    day = int(ele.text.strip())
                cols = row.find_all('td')
                cols = [ele.text.strip() for ele in cols] # a list of strings
                imo = 1
                for item in cols:
                    if len(item) == 0:
                        pass
                    else:
                        this_data = float(item.split()[0].replace(',',''))
                        # the split call is to remove trailing letters
                        # that occasionally appear after the data
                        d_dict[datetime(year, imo, day, 12, 0, 0)] = this_data
                    imo += 1
            rs.qt = pd.Series(d_dict)
            rs.flow_units = '$m^{3}s^{-1}$'
            #rs.fix_units() # not needed
            rs.qt = float(rs.ratio) * rs.qt
            rs.got_data = True
            rs.memo = 'success'
        else:
            rs.memo = 'That year was not available'
    except:
        rs.memo = 'problem parsing data from soup'
    rs.memo = (rs.memo + ' EC')

def fix_units(rs, qt):
    print(rs.flow_units)
    # fix units
    if rs.flow_units == 'kcfs':
        qt = qt*28.3168466
        rs.flow_units = 'm3/s'
    elif (rs.flow_units == 'cubic feet per second') or (rs.flow_units == 'ft3/s'):
        qt = qt*0.0283168466
        rs.flow_units = 'm3/s'
    else:
        pass
    return rs, qt