"""
Functions associated with NetCDF for tracker.
"""

from lo_tools import Lfun
import xarray as xr

# Info for NetCDF output, organized as {variable name: (long_name, units)}
name_unit_dict = {'lon':('Longitude','degrees'), 'lat':('Latitude','degrees'),
    'cs':('Fractional Z','Dimensionless'), 'ot':('Ocean Time',Lfun.roms_time_units),
    'z':('Z','m'), 'zeta':('Surface Z','m'), 'zbot':('Bottom Z','m'),
    'salt':('Salinity','Dimensionless'), 'temp':('Potential Temperature','Degrees C'),
    'oxygen':('Dissolved Oxygen', 'millimole_oxygen meter-3'),
    'u':('EW Velocity','meters s-1'), 'v':('NS Velocity','meters s-1'),
    'w':('Vertical Velocity','meters s-1'),
    'Uwind':('EW Wind Velocity','meters s-1'), 'Vwind':('NS Velocity','meters s-1'),
    'h':('Bottom Depth','m'),
    'hit_sidewall':('Hit Sidewall','1=hit'),
    'hit_bottom':('Hit Bottom','1=hit'),
    'hit_top':('Hit Top','1=hit'),
    'bad_pcs':('Bad Pcs','1=bad')}
    
def write_grid(g_infile, g_outfile):
    # write a file of grid info
    dsh = xr.open_dataset(g_infile)
    # lists of variables to process
    v_dict = dict()
    for vn in ['lon_rho', 'lat_rho', 'mask_rho', 'h']:
        v_dict[vn] = (('eta_rho', 'xi_rho'), dsh[vn].values)
    dsg = xr.Dataset(v_dict)
    dsg.to_netcdf(g_outfile)
    dsh.close()
    dsg.close()
    
def start_outfile(out_fn, P):
    v_dict = dict()
    for vn in P.keys():
        if vn == 'ot':
            v_dict[vn] = (('Time'), P[vn])
        else:
            v_dict[vn] = (('Time', 'Particle'), P[vn])
    ds = xr.Dataset(v_dict)
    for vn in P.keys():
        ds[vn].attrs['long_name'] = name_unit_dict[vn][0]
        ds[vn].attrs['units'] = name_unit_dict[vn][1]
    ds.to_netcdf(out_fn)
    ds.close()
    
def append_to_outfile(out_fn, P):
    """
    This works fine, but I worry about it when there are many days.  Will it slow down?
    
    It could be replaced by putting each daily output file in a temp dir, and then
    using ncrcat at the end.
    """
    v_dict = dict()
    for vn in P.keys():
        if vn == 'ot':
            v_dict[vn] = (('Time'), P[vn][1:])
        else:
            v_dict[vn] = (('Time', 'Particle'), P[vn][1:,:])
    ds1 = xr.Dataset(v_dict)
    ds0 = xr.open_dataset(out_fn, decode_times=False)
    ds2 = xr.concat((ds0,ds1), 'Time')
    ds0.close()
    ds1.close()
    out_fn.unlink(missing_ok=True)
    ds2.to_netcdf(out_fn)
    ds2.close()
