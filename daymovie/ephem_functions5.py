"""
Functions that make use of the ephem module.

"""
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import pytz
import ephem

AU2km = 149597871 # convert AU (astronomcal units) to km

def get_sunrise_sunset(dt0, dt1, city='Seattle'):
    """
    This returns lists of the datetimes sunrise and sunset
    at a given location (specified by the city argument).  The
    lists span the days given by the arguments dt0 and dt1 (datetimes).
    
    It calculates these in local time but then converts them to UTC.
    """
    if city == 'Seattle':
        city = 'Seattle'
        zone='US/Pacific'
    elif city == 'Westport':
        city = 'Westport'
        zone='US/Pacific'
    import ephem_functions as efun
    import pytz
    tz_utc, tz_local, obs = make_info(city=city, zone=zone)
    D_list = []
    D = dt0
    while D <= dt1:
        D_list.append(D)
        D += timedelta(days=1)
    Srise = []
    Sset = []
    for D in D_list:
        D_local = datetime(D.year, D.month, D.day, tzinfo=tz_local)
        S0, M0 = get_times(D_local, tz_utc, tz_local, obs)
        Srise.append(S0['rise'].astimezone(tz=pytz.timezone('UTC')))
        Sset.append(S0['set'].astimezone(tz=pytz.timezone('UTC')))
    return Srise, Sset

def make_info(city='Seattle', zone='US/Pacific'):
    # get timezone and observer info
    # use pytz.common_timezones to see a list of possibilities
    tz_local = pytz.timezone(zone)
    tz_utc = pytz.timezone('UTC')
    if city == 'Westport':
        # make obs by hand as a modification of Seattle
        obs = ephem.city('Seattle')
        obs.lon = '-124:06:36'
        obs.lat = '46:53:27'
    else:
        # or look if up in the ephem database (e.g. for Seattle)
        obs = ephem.city(city)
    return tz_utc, tz_local, obs

def get_sunmoon(dt_local, tz_utc, obs):
    # create sun and moon objects referenced to an observer
    # and a local datetime
    dt_utc = dt_local.astimezone(tz_utc)
    obs.date = ephem.Date(dt_utc)
    sun = ephem.Sun(obs)
    moon = ephem.Moon(obs)
    return sun, moon
    
def epht_to_dt(epht, tz=pytz.timezone('UTC')):
    # convert an ephem time to a datetime, including conversion
    # to a specific timezone if tz is provided (default = UTC)
    dt_naive = epht.datetime()
    dt_utc = pytz.utc.localize(dt_naive)
    dt = dt_utc.astimezone(tz)
    return dt
    
def get_times(dt_local, tz_utc, tz_local, obs, print_output=False):
    # gets rise, transit, and set times for sun and moon
    # in local timezone
    sun, moon = get_sunmoon(dt_local, tz_utc, obs)
    S = dict()
    M = dict()
    S['rise'] = epht_to_dt(obs.next_rising(sun), tz_local)
    S['transit'] = epht_to_dt(obs.next_transit(sun), tz_local)
    S['set'] = epht_to_dt(obs.next_setting(sun), tz_local)
    M['rise'] = epht_to_dt(obs.next_rising(moon), tz_local)
    M['transit'] = epht_to_dt(obs.next_transit(moon), tz_local)
    M['set'] = epht_to_dt(obs.next_setting(moon), tz_local)
    if print_output:
        fmt = '%14s %25s %s'
        tlist = ['rise', 'transit', 'set']
        print(fmt % ('Date', dt_local.ctime(), tz_local.zone))
        for tt in tlist:
            print(fmt % (('* sun '+tt), S[tt].ctime(), tz_local.zone))
        for tt in tlist:
            print(fmt % (('o moon '+tt), M[tt].ctime(), tz_local.zone))
    else:
        return S, M
        
def get_moon_orbit(dt0, dt1, daystep=1):
    # returns a pandas DataFrame with a time series of
    # distance, phase, and declination
    # from datetime dt0 to dt1 (assumed to be UTC)
    moon = ephem.Moon()
    dt_list = []
    dist_list = []
    phase_list = []
    dec_list = []
    dt = dt0
    while dt <= dt1:
        moon.compute(ephem.Date(dt))
        dt_list.append(dt)
        dist_list.append(moon.earth_distance * AU2km)
        # distance in AU = astronomical units = 149597871 km
        phase_list.append(moon.moon_phase)
        # phase as fraction of the surface illuminated
        dec_list.append(moon.dec * 180 / np.pi)
        # declination in radians
        dt = dt + timedelta(days=daystep)
    ddict = {'Distance (km)':dist_list,
            'Phase':phase_list,
            'Declination (deg)':dec_list}
    moon_orbit_df = pd.DataFrame(index=dt_list, data=ddict)
    return moon_orbit_df
    
def get_full_new(dt0, dt1):
    # get pandas DataFrames of full and new moons, including distance (km)
    # rounded to the nearest hour
    # from datetime dt0 to dt1 (assumed to be UTC)
    moon = ephem.Moon()
    fm_list = []
    dt = dt0
    dist_list = []
    while dt <= dt1:
        epht = ephem.Date(dt)
        epht = ephem.next_full_moon(epht)
        moon.compute(epht)
        dt = epht_to_dt(epht)
        if dt <= dt1:
            fm_list.append(dt)
            dist_list.append(moon.earth_distance * AU2km)
    full_df = pd.DataFrame(index=fm_list,
        data={'Distance (km)':dist_list})
    full_df.index = full_df.index.round('h')
    #
    nm_list = []
    dt = dt0
    dist_list = []
    while dt <= dt1:
        epht = ephem.Date(dt)
        epht = ephem.next_new_moon(epht)
        moon.compute(epht)
        dt = epht_to_dt(epht)
        if dt <= dt1:
            nm_list.append(dt)
            dist_list.append(moon.earth_distance * AU2km)
    new_df = pd.DataFrame(index=nm_list,
        data={'Distance (km)':dist_list})
    new_df.index = new_df.index.round('h')
    return full_df, new_df
    
if __name__ == '__main__':
    # example use of the functions
    tz_utc, tz_local, obs = make_info()
    dt_local = datetime(2016, 9, 19, tzinfo=tz_local)
    get_times(dt_local, tz_utc, tz_local, obs, print_output=True)
