"""
Code to read the ROMS text file varinfo.dat to construct the variable
dimensions and attributes used when writing to NetCDF.

"""

from lo_tools import Lfun

Ldir = Lfun.Lstart()

fn = Ldir['parent'] / 'LO_roms_source' / 'ROMS' / 'External' / 'varinfo.dat'

#vn_list = ['salt', 'temp', 'u', 'v', 'zeta', 'ubar', 'vbar']

dir_list = ['north', 'south', 'east', 'west']

vn_list = ['salt_' + item for item in dir_list]

# parse the file into a list of lists
ff = []
with open(fn, 'r') as f:
    for a in f:
        # make a list out of each line
        aa = a.replace('\n','').strip().split("'")
        aaa = [item for item in aa if len(item) > 0]
        ff.append(aaa)

v_dict = {}
for vn in vn_list:
    # make a dict for each variable in vn_list
    ii = 0
    for l in ff:
        if len(l) >= 2:
            # this line finds desired variables
            
            # this version works for state variables
            #if (l[0] == vn) and ('Input/Output' in l[1]):
                
            # this version works for boundary variables
            if (l[0] == vn):
                
                print('%d: %s' % (ii, l))
                # fill a dict for this variable based on the key below
                # copied from varinfo.dat
                d = {}
                d['long_name'] =    ff[ii+1][0]
                d['units'] =        ff[ii+2][0]
                d['field_type'] =   ff[ii+3][0]
                d['time_name'] =    ff[ii+4][0]
                d['index_name'] =   ff[ii+5][0]
                d['grid_type'] =    ff[ii+6][0]
                d['scale'] =        ff[ii+7][0]
                v_dict[vn] = d
        ii += 1
# !       Vinfo(1)      Field variable name.                                    !
# !       Vinfo(2)      Long-name attribute.                                    !
# !       Vinfo(3)      Units attribute.                                        !
# !       Vinfo(4)      Field type attribute.                                   !
# !       Vinfo(5)      Associated time variable name.                          !
# !       Vinfo(6)      Index variable name used in information arrays.         !
# !       Vinfo(7)      Staggered C-grid variable type:                         !
# !                       'nulvar' => non-grided variable.                      !
# !                       'p2dvar' => 2D PHI-variable.                          !
# !                       'r2dvar' => 2D RHO-variable.                          !
# !                       'u2dvar' => 2D U-variable.                            !
# !                       'v2dvar' => 2D V-variable.                            !
# !                       'p3dvar' => 3D PHI-variable.                          !
# !                       'r3dvar' => 3D RHO-variable.                          !
# !                       'u3dvar' => 3D U-variable.                            !
# !                       'v3dvar' => 3D V-variable.                            !
# !                       'w3dvar' => 3D W-variable.                            !
# !                       'b3dvar' => 3D BED-sediment.                          !
# !                       'l3dvar' => 3D spectral light variable.               !
# !                       'l4dvar' => 4D spectral light variable.               !
# !       scale         Scale to convert input data to model units.             !
