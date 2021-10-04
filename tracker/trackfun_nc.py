"""
Functions associated with NetCDF for tracker.
"""

import netCDF4 as nc4

# info for NetCDF output
name_unit_dict = {'lon':('Longitude','degrees'), 'lat':('Latitude','degrees'),
    'cs':('Fractional Z','Dimensionless'), 'ot':('Ocean Time','Seconds since 1/1/1970 UTC'),
    'z':('Z','m'), 'zeta':('Surface Z','m'), 'zbot':('Bottom Z','m'),
    'salt':('Salinity','Dimensionless'), 'temp':('Potential Temperature','Degrees C'),
    'u':('EW Velocity','meters s-1'), 'v':('NS Velocity','meters s-1'),
    'w':('Vertical Velocity','meters s-1'),
    'Uwind':('EW Wind Velocity','meters s-1'), 'Vwind':('NS Velocity','meters s-1'),
    'h':('Bottom Depth','m'),
    'hit_sidewall':('Hit Sidewall','1=hit'),
    'hit_bottom':('Hit Bottom','1=hit'),
    'hit_top':('Hit Top','1=hit'),
    'bad_pcs':('Bad Pcs','1=bad'),
}
    
def write_grid(g_infile, g_outfile):
    # write a file of grid info
    dsh = nc4.Dataset(g_infile)
    dsg = nc4.Dataset(g_outfile, 'w')
    # lists of variables to process
    dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
    # Copy dimensions
    for dname, the_dim in dsh.dimensions.items():
        if dname in dlist:
            dsg.createDimension(dname, len(the_dim))
    # Copy variables
    for vn in vn_list2:
        varin = dsh[vn]
        vv = dsg.createVariable(vn, varin.dtype, varin.dimensions)
        vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        vv[:] = dsh[vn][:]
    dsh.close()
    dsg.close()
    
def start_outfile(out_fn, P, NT_full, it0, it1):
    # This initializes for all anticipated time steps.
    NT, NP = P['lon'].shape
    ds = nc4.Dataset(out_fn, 'w')
    ds.createDimension('Time', NT_full)
    ds.createDimension('Particle', NP)
    # Copy variables
    for vn in P.keys():
        if vn == 'ot':
            vv = ds.createVariable(vn, float, ('Time'))
            vv[it0:it1] = P[vn][:]
        else:
            vv = ds.createVariable(vn, float, ('Time', 'Particle'))
            vv[it0:it1, :] = P[vn][:]
        vv.long_name = name_unit_dict[vn][0]
        vv.units = name_unit_dict[vn][1]
    ds.close()
    
def append_to_outfile(out_fn, P, it0, it1):
    # This is designed to add one day of data.
    ds = nc4.Dataset(out_fn, 'a')
    NTx, NPx = ds['lon'][:].shape
    for vn in P.keys():
        if vn == 'ot':
            ds[vn][it0:it1] = P[vn][:]
        else:
            ds[vn][it0:it1, :] = P[vn][:]
    ds.close()