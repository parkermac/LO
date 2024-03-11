"""
These are the functions we use to create our bio vars from salt.
Adapted from LO/forcing/ocn01/Ofun_bio.py.
"""

import numpy as np

def create_bio_var(salt, vn):
    # check salt limits
    if np.nanmax(salt) >= 35:
        print('Warning from create_bio_var(): salt > 35')
        print('Max salt value = %0.2f' % (np.nanmax(salt)))
        print('Action taken: limited salt in calculation to 34.99')
        salt[salt>=35] = 34.99
    
    if vn == 'NO3':
        # Salinity vs. NO3 [uM], Ryan McCabe 8/2015
        # NO3 = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 31.898)
        mm[ind] = 0
        bb[ind] = 0
        ind = ((salt >= 31.898) & (salt < 33.791))
        mm[ind] = 16.3958
        bb[ind] = -522.989
        ind = ((salt >= 33.791) & (salt < 34.202))
        mm[ind] = 29.6973
        bb[ind] = -972.4545
        ind = ((salt >= 34.202) & (salt < 34.482))
        mm[ind] = 8.0773
        bb[ind] = -233.0007
        ind = ((salt >= 34.482) & (salt < 35))
        mm[ind] = -28.6251
        bb[ind] = 1032.5686
        NO3 = mm*salt + bb
        # Set maximum NO3 to 45 microMolar (found at ~800m depth), based on
        # evidence from historical NO3 data in NODC World Ocean Database.
        NO3[NO3 > 45] = 45
        # Ensure that there are no negative values.
        NO3[NO3 < 0] = 0
        return NO3
    elif vn == 'NH4':
        NH4 = 0 * salt
        return NH4
    elif vn == 'oxygen':
        if np.nanmax(salt) > 36:
            print('Salt out of range for oxgen regression')
        # Salinity vs. oxygen [uM], Ryan McCabe 8/2015
        # oxygen = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 32.167)
        mm[ind] = 0
        bb[ind] = 300
        ind = ((salt >= 32.167) & (salt < 33.849))
        mm[ind] = -113.9481
        bb[ind] = 3965.3897
        ind = ((salt >= 33.849) & (salt < 34.131))
        mm[ind] = -278.3006
        bb[ind] = 9528.5742
        ind = ((salt >= 34.131) & (salt < 34.29))
        mm[ind] = -127.2707
        bb[ind] = 4373.7895
        ind = ((salt >= 34.29) & (salt < 34.478))
        mm[ind] = 34.7556
        bb[ind] = -1182.0779
        ind = ((salt >= 34.478) & (salt < 35))
        mm[ind] = 401.7916
        bb[ind] = -13836.8132
        oxygen = mm*salt + bb
        # Limit values.
        oxygen[oxygen > 450] = 450
        oxygen[oxygen < 0] = 0
        return oxygen
    elif vn == 'TIC':
        if np.nanmax(salt) > 36:
            print('Salt out of range for TIC regression')
        # Salinity vs. TIC [uM]
        # TIC = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 31.887)
        mm[ind] = 27.7967
        bb[ind] = 1112.2027
        ind = ((salt >= 31.887) & (salt < 33.926))
        mm[ind] = 147.002
        bb[ind] = -2688.8534
        ind = ((salt >= 33.926) & (salt < 34.197))
        mm[ind] = 352.9123
        bb[ind] = -9674.5448
        ind = ((salt >= 34.197) & (salt < 34.504))
        mm[ind] = 195.638
        bb[ind] = -4296.2223
        ind = ((salt >= 34.504) & (salt < 35))
        mm[ind] = -12.7457
        bb[ind] = 2893.77
        TIC = mm*salt + bb
        return TIC
    elif vn == 'alkalinity':
        mm = 0*salt
        bb = 0*salt
        if np.nanmax(salt) > 36:
            print('Salt out of range for alkalinity regression')
        # Salinity vs. alkalinity [uM]
        # alkalinity = mm*salt + bb
        ind = (salt < 31.477)
        mm[ind] = 37.0543
        bb[ind] = 1031.0726
        ind = ((salt >= 31.477) & (salt < 33.915))
        mm[ind] = 48.5821
        bb[ind] = 668.2143
        ind = ((salt >= 33.915) & (salt < 35))
        mm[ind] = 246.2214
        bb[ind] = -6034.6841
        alkalinity = mm*salt + bb
        return alkalinity
