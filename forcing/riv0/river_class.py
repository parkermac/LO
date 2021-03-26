"""
The river class code.
"""

import pandas as pd
import urllib
import xml.etree.ElementTree as ET
from datetime import datetime
#import numpy as np

class River:

    def __init__(self, riv_name, rs):
        # initialize all expected fields
        self.att_list = ['name', 'usgs_code',
            'nws_code', 'ec_code', 'scale_factor',
            'got_data', 'memo', 'flow_units']
        for att in self.att_list:
            setattr(self, att, '')
        self.name = riv_name
        self.rs = rs
        # default values
        self.Q = []
        self.T = []
        self.qt = pd.Series(self.Q, index=self.T)
        try:
            self.scale_factor = self.rs['ratio']
            self.usgs_code = self.rs['usgs']
            self.nws_code = self.rs['nws']
            self.ec_code = self.rs['ec']
        except KeyError:
            # needed for analytical rivers
            pass
        self.got_data = False
        self.memo = 'no message'

    def print_info(self):
        print(50*'-')
        for att in self.att_list:
            print(att + ' = ' + str(getattr(self,att)))
        print(50*'-')
        
    def get_usgs_data(self, days):
        # this deals with rivers that need gauge combinations
        if self.name == 'skokomish':
            print('+++ combining to form Skokomish River +++')
            self.usgs_code = 12060500
            self.scale_factor = 1.4417
            self.get_usgs_data_sub(days)
            qt1 = self.qt.copy()
            self.usgs_code = 12059500
            self.scale_factor = 1.0
            self.get_usgs_data_sub(days)
            qt2 = self.qt.copy()
            self.qt = qt1 + qt2
        elif self.name == 'hamma':
            print('+++ combining to form Hamma Hamma River +++')
            self.usgs_code = 12060500
            self.scale_factor = 1.4417
            self.get_usgs_data_sub(days)
            qt1 = self.qt.copy()
            self.usgs_code = 12059500
            self.scale_factor = 1.0
            self.get_usgs_data_sub(days)
            qt2 = self.qt.copy()
            self.qt = 0.4125 * (qt1 + qt2)
        else:
            self.get_usgs_data_sub(days)

    def get_usgs_data_sub(self, days):
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
        gage = int(self.usgs_code)
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
            for e0 in root.findall(".//"):
                if e0.tag == aa+'value':
                    self.Q.append(float(e0.text))
                    self.T.append(pd.to_datetime(e0.get('dateTime')))
                if e0.tag == aa+'unitCode' and flag:
                    self.flow_units = e0.text
                    flag = False
            self.qt = pd.Series(self.Q, index=self.T)
            self.fix_units()
            self.qt = float(self.scale_factor) * self.qt
            # Note: this resampling of daily data just moves the timestamp to noon
            # of the day it started at.  Data is unchanged.
            self.qt = self.qt.resample('D', label='right', loffset='-12h').mean()
            self.got_data = True
            self.memo = 'success'
        except:
            self.memo = 'problem downloading or parsing data from XML'
        self.memo = (self.memo + ' USGS')

    def get_nws_data(self):
        # This gets NWS forecast data.
        url_str = ('http://www.nwrfc.noaa.gov/xml/xml.cgi?id=' +
                   self.nws_code +
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
                                    self.Q.append(float(ee.text))
                                    if flag:
                                        self.flow_units = ee.get('units')
                                        flag = False
                                if ee.tag == aa+'dataDateTime':
                                    self.T.append(pd.to_datetime(ee.text))
            self.qt = pd.Series(self.Q, index=self.T)
            self.fix_units()
            self.qt = float(self.scale_factor) * self.qt
            self.qt = self.qt.resample('D', label='right', loffset='-12h').mean()
            self.got_data = True
            self.memo = 'success'
        except:
            self.memo = 'problem downloading or parsing data from XML'
        self.memo = (self.memo + ' NWS')

    def get_ec_data(self, days, timeout=10):
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
                'stn': self.ec_code,
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
            self.qt = pd.Series(d_dict)
            self.flow_units = '$m^{3}s^{-1}$'
            #self.fix_units() # not needed
            self.qt = float(self.scale_factor) * self.qt
            self.qt = self.qt.resample('D', label='right', loffset='-12h').mean()
            # NEW 2019.03.20 to deal with the problem that when you request date from
            # before 18 months ago it gives the most recent data instead.
            dt0_actual = self.qt.index[0]
            dt0_requested = days[0]
            import numpy as np
            if np.abs((dt0_actual - dt0_requested).days) >= 1:
                memo = 'That date range was not available'
                qt = ''
            else:
                self.got_data = True
                self.memo = 'success'
        except:
            self.memo = 'problem parsing data from soup'
        self.memo = (self.memo + ' EC')

    def get_ec_data_historical(self, year):
        # gets Environment Canada data, using code cribbed from:
        #https://bitbucket.org/douglatornell/ecget/src/
        # NOTE: this will get data up through the end of 2016.
        import requests
        import bs4
        try:
            params = {
                'mode': 'Table',
                'type': 'h2oArc',
                'stn': self.ec_code,
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
                self.qt = pd.Series(d_dict)
                self.flow_units = '$m^{3}s^{-1}$'
                #self.fix_units() # not needed
                self.qt = float(self.scale_factor) * self.qt
                self.got_data = True
                self.memo = 'success'
            else:
                self.memo = 'That year was not available'
        except:
            self.memo = 'problem parsing data from soup'
        self.memo = (self.memo + ' EC')

    def fix_units(self):
        # fix units
        if self.flow_units == 'kcfs':
            self.qt = self.qt*28.3168466
            self.flow_units = '$m^{3}s^{-1}$'
        elif (self.flow_units == 'cubic feet per second'
            or self.flow_units == 'ft3/s'):
            self.qt = self.qt*0.0283168466
            self.flow_units = '$m^{3}s^{-1}$'
        else:
            pass


