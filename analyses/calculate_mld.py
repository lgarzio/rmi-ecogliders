#!/usr/bin/env python

"""
Author: Lori Garzio on 1/11/2022
Last modified: 6/24/2025
Modified from code from Sam Coakley following theory from Carvalho et al 2016 https://doi.org/10.1002/2016GL071205
Calculate Mixed Layer Depth for glider profiles using density and pressure, then add the MLD variable to the .nc file.
The dataset provided must have 'time' or 'profile_time' as the only coordinate in order to convert the dataset to a
dataframe properly
"""

import os
import datetime as dt
import xarray as xr
import numpy as np
import pandas as pd
import gsw
import matplotlib.pyplot as plt
import functions.common as cf
import functions.mixed_layer_depth as mldfunc
plt.rcParams.update({'font.size': 14})
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


def main(fname, timevar, mldvar, zvar, qithresh):
    mld = np.array([], dtype='float32')
    max_n2 = np.array([], dtype='float32')

    mld_summary = 0
    mld_summary_nan = 0
    mld_summary_lessthan_qi = 0

    savefile = f'{fname.split(".nc")[0]}_mld.nc'

    ds = xr.open_dataset(fname)
    ds = ds.sortby(ds.time)
    try:
        deploy = ds.attrs['deployment']
    except KeyError:
        f = fname.split('/')[-1]
        deploy = f'{f.split("-")[0]}-{f.split("-")[1]}'  # get the deployment from the filename

    df = ds.to_dataframe()

    # remove data that's collected at the surface (< 1 dbar/m)
    df.loc[df[zvar] < 1, [zvar, mldvar]] = np.nan
    grouped = df.groupby(timevar, dropna=False)

    for group in grouped:
        # Create temporary dataframe to interpolate to dz m depths
        ll = len(group[1])
        kwargs = {'depth_var': zvar}
        temp_df1 = group[1][[mldvar, zvar]].dropna(how='all')
        if len(temp_df1) == 0:
            mldx = np.repeat(np.nan, ll)
            max_n2x = np.repeat(np.nan, ll)
        else:
            temp_df = cf.depth_bin(temp_df1, **kwargs)
            temp_df.dropna(subset=[mldvar], inplace=True)
            temp_df.index.name = f'{zvar}_bins'
            temp_df.reset_index(inplace=True)
            if len(temp_df) == 0:
                mldx = np.repeat(np.nan, ll)
                max_n2x = np.repeat(np.nan, ll)
                qi = 'MLD not calculated'
            else:
                # calculate profile's pressure range
                pressure_range = (np.nanmax(temp_df[zvar]) - np.nanmin(temp_df[zvar]))

                if pressure_range < 5:
                    # if the profile spans <5 dbar, don't calculate MLD
                    mldx = np.repeat(np.nan, ll)
                    max_n2x = np.repeat(np.nan, ll)
                    qi = 'MLD not calculated'

                else:
                    #kwargs = {'zvar': zvar}
                    kwargs = {'zvar': zvar, 'qi_threshold': qithresh}
                    mldx, N2, qi = mldfunc.profile_mld(temp_df, **kwargs)
                    mldx = np.repeat(mldx, ll)
                    max_n2x = np.repeat(N2, ll)

        mld = np.append(mld, mldx)
        max_n2 = np.append(max_n2, max_n2x)

        # summarize the MLD results
        if np.sum(~np.isnan(mldx)) > 0:
            mld_summary += 1
        elif isinstance(qi, float):
            mld_summary_lessthan_qi += 1
        else:
            mld_summary_nan += 1

    # add mld to the dataset
    mld_min = np.round(np.nanmin(mld), 4)
    mld_max = np.round(np.nanmax(mld), 4)
    attrs = {
        'actual_range': np.array([mld_min, mld_max]),
        'ancillary_variables': [mldvar, zvar],
        'observation_type': 'calculated',
        'units': ds[zvar].units,
        'comment': 'Mixed Layer Depth for each profile, calculated as the depth of max Brunt‐Vaisala frequency squared (N**2) '
                   'from Carvalho et al 2016 (https://doi.org/10.1002/2016GL071205)',
        'long_name': 'Mixed Layer Depth'
        }
    da = xr.DataArray(mld, coords=ds[mldvar].coords, dims=ds[mldvar].dims,
                      name='mld', attrs=attrs)
    ds['mld'] = da

    # add maximum buoyancy frequency N2 (measure of stratification strength) to the dataset
    n2_min = np.round(np.nanmin(max_n2), 6)
    n2_max = np.round(np.nanmax(max_n2), 6)
    attrs = {
        'actual_range': np.array([n2_min, n2_max]),
        'ancillary_variables': [mldvar, zvar],
        'observation_type': 'calculated',
        'units': 's-2',
        'comment': 'Maximum Brunt‐Vaisala frequency squared (N**2) for each profile used to calculate Mixed Layer Depth '
                   'from Carvalho et al 2016 (https://doi.org/10.1002/2016GL071205). This can be used as a measurement '
                   'for stratification strength',
        'long_name': 'Maximum Buoyancy Frequency'
    }
    da = xr.DataArray(max_n2, coords=ds[mldvar].coords, dims=ds[mldvar].dims,
                      name='max_n2', attrs=attrs)
    ds['max_n2'] = da

    ds.to_netcdf(savefile)

    print(deploy)
    print(f'Total profiles: {len(grouped)}')
    print(f'MLD calculated for {mld_summary} profiles')
    print(f'MLD not calculated for {mld_summary_lessthan_qi} profiles because QI < {qithresh}')
    print(f'MLD not calculated for {mld_summary_nan} profiles for other reasons (e.g. pressure range < 5 dbar, no data)')


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru39-20241021T1717/ncei_pH/ru39-20241021T1717-delayed.nc'
    time_variable = 'profile_time'  # time variable on which groups are generated (e.g. profile_time)
    mldvar = 'density'  # variable used to calculate MLD
    zvar = 'depth'  # pressure variable
    qithreshold = 0.8  # 0.5 0.8
    main(ncfile, time_variable, mldvar, zvar, qithreshold)
