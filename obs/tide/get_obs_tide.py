"""
Code to automate getting year-long tide height records from
a series of NOAA and DFO sites around the Salish Sea and NE Pacific
coast.

Also computes harmonics using utide.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun

import pickle

from importlib import reload
import obsfun as ofn
reload(ofn)

home = os.environ.get('HOME')
dir00 = home + '/Documents/'
dir0 = dir00 + 'ptools_output/tide/'

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()

# extract and save data
year  = 2017
outdir = dir0 + 'obs_data/'
Lfun.make_dir(outdir)

testing = False
if testing == True:
    noaa_sn_dict = {
        'Charleston': 9432780}
    dfo_sn_dict = {
        'Point Atkinson': 7795}

for name in noaa_sn_dict.keys():
    sn = noaa_sn_dict[name]
    fn = outdir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = outdir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = outdir + 'h_' + str(sn) + '_' + str(year) + '.p'
    print(name)
    df, m_dict = ofn.get_noaa_tide(sn, year)
    h = ofn.get_harmonics(df, float(m_dict['lat']))
    df.to_pickle(fn)
    Lfun.dict_to_csv(m_dict, mfn)
    pickle.dump(h, open(hfn, 'wb'))

for name in dfo_sn_dict.keys():
    sn = dfo_sn_dict[name]
    fn = outdir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = outdir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = outdir + 'h_' + str(sn) + '_' + str(year) + '.p'
    print(name)
    df, m_dict = ofn.get_dfo_tide(sn, year)
    h = ofn.get_harmonics(df, float(m_dict['lat']))
    df.to_pickle(fn)
    Lfun.dict_to_csv(m_dict, mfn)
    pickle.dump(h, open(hfn, 'wb'))
    