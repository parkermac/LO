"""
Functions to work with observational tide data.
"""

import requests
import xml.etree.ElementTree as ET
import pandas as pd
from datetime import datetime
import pytz
# import utide
from matplotlib.dates import date2num

def get_AG(hn, Hobs, Hmod):
    #convenience function for loading constituent info
    ho = Hobs
    hm = Hmod
    # we use the "[0]" because these are arrays and we want floats
    Ao = ho.A[ho.name==hn][0]
    Am = hm.A[hm.name==hn][0]
    Go = ho.g[ho.name==hn][0]
    Gm = hm.g[hm.name==hn][0]
    Fo = 24*ho.aux.frq[ho.name==hn][0] # cycles per day
    Fm = 24*hm.aux.frq[hm.name==hn][0]
    # fix when phase difference straddles 360
    if (Gm - Go) > 180:
        Gm = Gm - 360
    elif (Gm - Go) < -180:
        Gm = Gm + 360
    else:
        pass
    return Ao, Am, Go, Gm, Fo, Fm

# list of frequencies to consider.  Sometimes we want to limit this
# because for shorter records utide can't separate nearby peaks
#hn_list = ['M2','S2','N2','O1','P1','K1']
hn_list = ['M2','S2','N2','O1','K1']

def get_sn_dicts():
    noaa_sn_dict = {
        'Charleston': 9432780,
        'South Beach': 9435380,
        'Garibaldi': 9437540,
        'Toke Point': 9440910,
        'Westport': 9441102,
        'La Push': 9442396,
        'Neah Bay': 9443090,
        'Port Angeles': 9444090,
        'Friday Harbor': 9449880,
        'Cherry Point': 9449424,
        'Port Townsend': 9444900,
        'Seattle': 9447130,
        'Tacoma': 9446484}
    dfo_sn_dict = {
        'Point Atkinson': 7795,
        'Vancouver': 7735,
        'Patricia Bay': 7277,
        'Victoria Harbour': 7120,
        'Bamfield': 8545,
        'Tofino': 8615,
        #'Winter Harbour': 8735,
        #'Port Hardy': 8408,
        'Campbell River': 8074,
        'New Westminster': 7654}
    sn_dict = {}
    sn_dict.update(noaa_sn_dict)
    sn_dict.update(dfo_sn_dict)
    #
    return noaa_sn_dict, dfo_sn_dict, sn_dict

def get_noaa_tide(sn, year):
    # inputs can be ints or strings
    sn = str(sn) # station number
    year_str = str(year)
    # next_year_str = str(year+1)
    a = ('https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?'
        + 'begin_date=' + year_str + '0101'
        + '&end_date=' + year_str + '1231'
        + '&station=' + sn
        + '&product=hourly_height'
        + '&datum=MLLW&units=metric&time_zone=gmt'
        + '&application=University_of_Washington'
        + '&format=xml')
    b = requests.get(a)
    root = ET.fromstring(b.text)
    # metadata
    m_dict = dict()
    m = root.find('metadata')
    for key in m.keys():
        m_dict[key] = m.attrib[key]
    # data
    t_list = []
    eta_list = []
    for e0 in root.findall('observations'):
        for e in e0.findall('hr'):
            t_list.append(e.attrib['t'])
            eta_list.append(float(e.attrib['v']))
    dti = pd.to_datetime(t_list)
    dti = dti.tz_localize('UTC')
    df = pd.DataFrame(data={'eta':eta_list}, index = dti)
    df.index.name = 'Date'
    return df, m_dict
    
def get_dfo_tide(sn, year):
    sn = str(sn) # station number
    year = str(year)
    start_date = '31-DEC-' + str(int(year)-1)
    end_date = '01-JAN-' + str(int(year)+1)
    outfile = '../../ptools_data/tide/dfo_scratch_'+sn+'_'+year+'.csv'
    # Form urls and html information
    base_url = 'http://www.meds-sdmm.dfo-mpo.gc.ca/isdm-gdsi/twl-mne/inventory-inventaire/'
    form_handler = ('data-donnees-eng.asp?user=isdm-gdsi&region=PAC&tst=1&no='
        + sn)
    sitedata = {'start_period': start_date,
        'end_period': end_date,
        'resolution': 'h',
        'time_zone': 'u'}
    data_provider = (
        'download-telecharger.asp'
        '?File=E:%5Ciusr_tmpfiles%5CTWL%5C'
        + sn + '-'+start_date + '_slev.csv'
        '&Name=' + sn + '-'+start_date+'_slev.csv')
    # Go get the data from the DFO site
    with requests.Session() as s:
        s.post(base_url + form_handler, data=sitedata)
        r = s.get(base_url + data_provider)
    # Write the data to a text file
    with open(outfile, 'w') as f:
        f.write(r.text)
    # and parse the text file
    df = read_dfo_tide(outfile, year)
    m_dict = read_dfo_info(outfile)
    return df, m_dict
    
def read_dfo_tide(fn, year):
    df = pd.read_csv(fn, skiprows=8, header=None)
    df.columns = ['Date','eta','junk']
    df = df.set_index('Date')
    df = df.drop(['junk'], axis=1)
    df.index = pd.to_datetime(df.index)
    df = df.tz_localize('UTC')
    year = int(year)
    dt0 = datetime(year,1,1,tzinfo=pytz.timezone('UTC'))
    dt1 = datetime(year,12,31,23,tzinfo=pytz.timezone('UTC'))
    df = df[dt0:dt1]
    return df

def read_dfo_info(fn):
    df = pd.read_csv(fn, nrows=6, index_col=0, header=None)
    df.index.name = 'Item'
    df.columns = ['Value', 'junk']
    df = df.drop(['junk'], axis=1)
    mm_dict = df['Value'].to_dict()
    # translate keys
    in_list = ['Station_Name', 'Station_Number',
        'Latitude_Decimal_Degrees', 'Longitude_Decimal_Degrees']
    out_list = ['name', 'id', 'lat', 'lon']
    name_dict = dict(zip(in_list, out_list))
    m_dict = dict()
    for key in in_list:
        m_dict[name_dict[key]] = mm_dict[key]
    m_dict['lon'] = '-' + m_dict['lon']
    return m_dict
    
def get_harmonics(df, lat):
    t = date2num(df.index.to_pydatetime())
    # NOTE 1/23/2025 t should just be a DatetimeIndex
    z = df['eta'].to_numpy()
    h = utide.solve(t, z, v=None,
                 lat=lat,
                 nodal=False,
                 trend=False,
                 method='ols',
                 conf_int='linear',
                 Rayleigh_min=0.95)
    # h.aux.freq has units cyles/hour
    # so for f = h.aux.frq[h.name == 'M2'][0] we get
    # 1/f = 12.420601202671868 (hours per cycle)
    # h.A is amplitude (m), h.g is phase (degrees)
    return h
             
if __name__ == '__main__':
    # examples of uses of the functions
    if True:
        # noaa
        sn = 9447130
        year = 2013
        df, m_dict = get_noaa_tide(sn, year)
        lat=float(m_dict['lat'])
        h = get_harmonics(df, lat)
    else:
        # dfo
        sn = 7795
        year = 2013
        df, m_dict = get_dfo_tide(sn, year)
        #
        lat=float(m_dict['lat'])
        h = get_harmonics(df, lat)
    if True:
        import matplotlib.pyplot as plt
        plt.close('all')
        df.plot(title=m_dict['name'])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for ii in range(len(h.A)):
            # plot amplitude vs. frequency
            ax.text(h.aux.frq[ii], h.A[ii], h.name[ii])
        ax.set_title(m_dict['name'])
        plt.show()
        
