#!/usr/bin/env python

"""
Author: Lori Garzio on 4/3/2025
Last modified: 4/3/2025
Modified from code from Julia Engdahl
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


def main(fname, timevar):
    mld_upper = np.array([], dtype='float32')
    mld_lower = np.array([], dtype='float32')
    qc_upper = np.array([], dtype='float32')
    qc_lower = np.array([], dtype='float32')

    savefile = f'{fname.split("_mld.nc")[0]}_mld2way.nc'

    ds = xr.open_dataset(fname)
    ds = ds.sortby(ds.time)
    try:
        deploy = ds.attrs['deployment']
    except KeyError:
        f = fname.split('/')[-1]
        deploy = f'{f.split("-")[0]}-{f.split("-")[1]}'  # get the deployment from the filename

    df = ds.to_dataframe()

    mldvars = ['temperature', 'salinity']
    select_vars = [timevar, 'profile_lat', 'profile_lon', 'depth_interpolated', 'temperature', 'salinity']

    grouped = df.groupby(timevar, dropna=False)

    for group in grouped:
        # Create temporary dataframe to interpolate to dz m depths
        ll = len(group[1])
        temp_df = group[1][select_vars].dropna(how='all', subset=mldvars)
        if len(temp_df) == 0:
            mldx_upper = np.repeat(np.nan, ll)
            mldx_lower = np.repeat(np.nan, ll)
        else:
            # calculate profile's depth range
            depth_range = np.nanmax(temp_df['depth_interpolated']) - np.nanmin(temp_df['depth_interpolated'])

            if depth_range < 10:
                # if the profile spans <10 m, don't calculate MLD, qc=4
                mldx_upper = np.repeat(np.nan, ll)
                mldx_lower = np.repeat(np.nan, ll)
                qcx_upper = np.repeat(4, ll)
                qcx_lower = np.repeat(4, ll)

            else:
                # sort the dataframe by depth
                temp_df = temp_df.sort_values(by='depth_interpolated')

                mldU, mldL, qcU, qcL = mldfunc.profile_mld_2way(temp_df)
                mldx_upper = np.repeat(mldU, ll)
                mldx_lower = np.repeat(mldL, ll)
                qcx_upper = np.repeat(qcU, ll)
                qcx_lower = np.repeat(qcL, ll)

                if mldU < 3:
                    print('mldU < 3')
                    #mldU, mldL, qcU, qcL = mldfunc.profile_mld_2way(temp_df)
                if mldL < 3:
                    print('mldL < 3')
                    #mldU, mldL, qcU, qcL = mldfunc.profile_mld_2way(temp_df)

        mld_upper = np.append(mld_upper, mldx_upper)
        mld_lower = np.append(mld_lower, mldx_lower)
        qc_upper = np.append(qc_upper, qcx_upper)
        qc_lower = np.append(qc_lower, qcx_lower)

        # # summarize the MLD results
        # if np.sum(~np.isnan(mldx)) > 0:
        #     mld_summary += 1
        # elif isinstance(qi, float):
        #     mld_summary_lessthan_qi += 1
        # else:
        #     mld_summary_nan += 1

    # add mld upper to the dataset
    mldu_min = np.round(np.nanmin(mld_upper), 4)
    mldu_max = np.round(np.nanmax(mld_upper), 4)
    attrs = {
        'actual_range': np.array([mldu_min, mldu_max]),
        'observation_type': 'calculated',
        'units': 'm',
        'comment': 'Upper boundary of the transition zone for mixed layer depth. Based off of de Boyer Montegut (2007) and Rudzin (2017)',
        'long_name': 'Mixed Layer Depth Upper',
        }
    da = xr.DataArray(mld_upper, coords=ds['salinity'].coords, dims=ds['salinity'].dims,
                      name='mld_upper', attrs=attrs)
    ds['mld_upper'] = da

    # add mld lower to the dataset
    mldl_min = np.round(np.nanmin(mld_lower), 4)
    mldl_max = np.round(np.nanmax(mld_lower), 4)
    attrs = {
        'actual_range': np.array([mldl_min, mldl_max]),
        'observation_type': 'calculated',
        'units': 'm',
        'comment': 'Lower boundary of the transition zone for mixed layer depth. Based off of de Boyer Montegut (2007) and Rudzin (2017)',
        'long_name': 'Mixed Layer Depth Lower',
        }
    da = xr.DataArray(mld_lower, coords=ds['salinity'].coords, dims=ds['salinity'].dims,
                      name='mld_upper', attrs=attrs)
    ds['mld_lower'] = da

    ds.to_netcdf(savefile)

    print(deploy)
    # print(f'Total profiles: {len(grouped)}')
    # print(f'MLD calculated for {mld_summary} profiles')
    # print(f'MLD not calculated for {mld_summary_lessthan_qi} profiles because QI < {qithresh}')
    # print(f'MLD not calculated for {mld_summary_nan} profiles for other reasons (e.g. pressure range < 5 dbar, no data)')


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2023/ru39-20230817T1520/ncei_pH/ru39-20230817T1520-delayed_mld.nc'
    time_variable = 'profile_time'  # time variable on which groups are generated (e.g. profile_time)
    main(ncfile, time_variable)
