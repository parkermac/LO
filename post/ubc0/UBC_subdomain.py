# Python script for extracting a subset of Live Ocean results
# ---------------------------------------------------------------------------
# Usage
#
# From command line:
# python UBC_subdomain.py filename1 filename2 filename3 ...
#
# Inside python:
# import UBC_subdomain
# UBC_subdomain.get_UBC_subdomain([filename1, filename2, filename3,...])
# ----------------------------------------------------------------------------
# Description
#
# Creates new netCDF files that store a subdomain and subset of variables
# for UBC users (S. Allen's group).
# New file names are derived from the original file name with '_UBC' suffix.
# For example: ocean_his_0002.nc becomes ocean_his_0002_UBC.nc
# ----------------------------------------------------------------------------
# Written by
# Nancy Soontiens 2016
# nsoontie@eos.ubc.ca
#
# Modified 2018.11.08 by P. MacCready to get the indices from the calling function
# make_forcing_main.py instead of hard-coding them.  This means the the lat-lon
# range to extract is set in the calling function.  I also changed the format of
# the netcdf ourput to be NETCDF4_CLASSIC so that is works with the call to
# MFDataset used by the low pass code in the calling function.

import sys
import netCDF4 as nc

ncformat = 'NETCDF3_64BIT_OFFSET' # PM Edit 2019.05.15


# XBS = [55, 80]  # x-limits
# YBS = [295, 325]  # y-limits
# Variables to copy
VAR_LIST = ['salt', 'temp', 'h', 'lon_rho', 'lat_rho', 'mask_rho', 'pn', 'pm',
            's_rho', 'hc', 'Cs_r', 'Vtransform', 'zeta', 'ocean_time',
            'lon_u', 'lat_u', 'mask_u', 'u',
            'lon_v', 'lat_v', 'mask_v', 'v',
            'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
            'oxygen', 'TIC', 'alkalinity', 'CaCO3', 'rho']
# Dimensions to copy
# DIM_LIST = ['xi_rho', 'eta_rho', 'N', 's_rho', 'ocean_time',
#             'xi_u', 'eta_u', 'xi_v', 'eta_v']
DIM_LIST = ['xi_rho', 'eta_rho', 's_rho', 'ocean_time',
            'xi_u', 'eta_u', 'xi_v', 'eta_v']


def get_UBC_subdomain(f_list, out_dir, XBS, YBS):
    """Create subdomain files for all netCDF files in f_list """
    for fname in f_list:
        fnew = out_dir / fname.name.replace('.nc', '_UBC.nc')
        #print(fnew)
        #sys.stdout.flush()
        #fnew = '{}_UBC.nc'.format(fname.split('.nc', 1)[0])
        # with nc.Dataset(fname) as G, nc.Dataset(fnew, 'w', format='NETCDF4_CLASSIC') as Gnew:
        #     _copy_netCDF_subdomain(G, Gnew, XBS, YBS, VAR_LIST, DIM_LIST)
        with nc.Dataset(fname) as G, nc.Dataset(fnew, 'w', format=ncformat) as Gnew:
            _copy_netCDF_subdomain(G, Gnew, XBS, YBS, VAR_LIST, DIM_LIST)


def _copy_netCDF_subdomain(oldfile, newfile, xbounds, ybounds,
                           var_list, dim_list):
    """Copy variables in var_list in subdomain [xbounds, ybounds] from
       oldfile to newfile. Also copies dimensions in dim_list and
       all global attributes.
    """
    _copy_dimensions(oldfile, newfile, dim_list, xbounds, ybounds)
    _copy_variables(oldfile, newfile, var_list, xbounds, ybounds)
    # copy global attributes
    newfile.setncatts(
        {att: oldfile.getncattr(att) for att in oldfile.ncattrs()}
        )


def _copy_dimensions(oldfile, newfile, dim_list, xbounds, ybounds):
    """ Copy the dimensions in dims_list from oldfile to newfile.
        Dimensions of eta_rho, xi_rho are determined by limits of
        ybounds, xbounds. eta_v and xi_u have one extra because of staggering.
    """
    dim_size_dict = {
        'eta_rho': ybounds[1]-ybounds[0]+1,
        'eta_u':  ybounds[1]-ybounds[0]+1,
        'xi_rho': xbounds[1]-xbounds[0]+1,
        'xi_v': xbounds[1]-xbounds[0]+1,
        'xi_u': xbounds[1]-xbounds[0],
        'eta_v': ybounds[1]-ybounds[0],
        'ocean_time': 0,
    }
    for dimname in dim_list:
        dim = oldfile.dimensions[dimname]
        newfile.createDimension(
            dimname, size=dim_size_dict.get(dimname, dim.__len__())
        )


def _copy_variables(oldfile, newfile, var_list, xbounds, ybounds):
    """Copy variables in var_list from oldfile to newfile for subdomain
        [xbounds, ybounds]"""
    varnames_in_file = list(oldfile.variables.keys())
    for varname in var_list:
        if varname in varnames_in_file:
            var = oldfile.variables[varname]
            dims = var.dimensions
            newvar = newfile.createVariable(varname, var.datatype, dims)
            # copy variable attributes
            newvar.setncatts({att: var.getncattr(att)
                             for att in var.ncattrs()})
            # fill data
            if 'eta_rho' in dims or 'xi_rho' in dims:
                newvar[:] = var[...,
                                ybounds[0]:ybounds[1]+1,
                                xbounds[0]:xbounds[1]+1]
            elif 'eta_u' in dims or 'xi_u' in dims:
                newvar[:] = var[...,
                                ybounds[0]:ybounds[1]+1,
                                xbounds[0]:xbounds[1]]
            elif 'eta_v' in dims or 'xi_v' in dims:
                newvar[:] = var[...,
                                ybounds[0]:ybounds[1],
                                xbounds[0]:xbounds[1]+1]
            else:
                newvar[:] = var[:]


if __name__ == '__main__':
    get_UBC_subdomain(sys.argv[1:])
