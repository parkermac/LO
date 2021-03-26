# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 08:52:07 2016

@author: PM5
"""

"""
The river class code.  This version gets Temperature instead of flow.

We still use Q and qt for the fields, but Q or q refer to temperature in
this case. (yuck!)

"""

import pandas as pd
import urllib
import xml.etree.ElementTree as ET

class River:

    def __init__(self, riv_name, rs):
        # initialize all expected fields
        self.att_list = ['name', 'usgs_code', 'ec_code',
                         'got_data', 'memo', 'temp_units']
        for att in self.att_list:
            setattr(self, att, '')
        self.name = riv_name
        self.rs = rs
        # default values
        self.Q = []
        self.T = []
        self.qt = pd.Series(self.Q, index=self.T)
        try:
            self.usgs_code = self.rs['usgs']
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
            + time_str + '&parameterCd=00010')
        # Get the XML.
        try: # an extra try/except loop
            file = urllib.request.urlopen(url_str, timeout=10)
        except:
            self.memo = 'problem opening URL'
        try:
            tree = ET.parse(file)
            root = tree.getroot()
        except:
            self.memo = 'problem parsing XML'
        # Extract the data from the XML.
        try:
            flag = True
            # aa = '{http://www.cuahsi.org/waterML/1.1/}'
            rt = root.tag
            aa = rt[rt.find('{'): rt.find('}') + 1]
            for e0 in root.findall(".//"):
                if e0.tag == aa+'value':
                    self.Q.append(float(e0.text))
                    self.T.append(pd.to_datetime(e0.get('dateTime')))
                if e0.tag == aa+'unitCode' and flag:
                    self.temp_units = e0.text
                    flag = False
            self.qt = pd.Series(self.Q, index=self.T)
            self.qt = self.qt.resample('D', label='right', loffset='-12h').mean()
            self.got_data = True
            self.memo = 'success'
        except:
            self.memo = 'problem getting data from XML'
        self.memo = (self.memo + ' USGS')

    def get_ec_data(self, days, timeout=10):
        # NOTE: this will get data up through today, but can only go back
        # 18 months into the past.

        import requests
        import bs4
        try:
            PARAM_IDS = {'discharge': 47, 'temperature': 5}
            params = {
                'mode': 'Table',
                'type': 'realTime',
                'prm1': PARAM_IDS['temperature'],
                'prm2': -1,
                'stn': self.ec_code,
                'startDate': days[0].strftime('%Y-%m-%d'),
                'endDate': days[1].strftime('%Y-%m-%d'),
            }
            DATA_URL = 'http://wateroffice.ec.gc.ca/report/report_e.html'
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
                d_dict[pd.to_datetime(item[0])] = float(item[1])
            self.qt = pd.Series(d_dict)
            self.temp_units = 'deg C'
            self.qt = self.qt.resample('D', label='right', loffset='-12h').mean()
            self.got_data = True
            self.memo = 'success'
        except:
            self.memo = 'problem parsing data from soup'
        self.memo = (self.memo + ' EC')


