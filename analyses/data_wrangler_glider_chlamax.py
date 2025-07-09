#!/usr/bin/env python

"""
Author: Lori Garzio on 7/1/2025
Last modified: 7/8/2025
Glider data wrangler, export as NetCDF. The netcdf files have already been QC'd and processed, so this script
is for synthesizing data at the chl maximum.
For each profile: if a mixed layer depth (N2 buoyancy frequency version) is calculated, then
1. Bin data into 1m depth bins
1. Find the chla max for each binned profile
2. Find the depth of the chla max
3. Grab the pH and DO data in the +/- 3m depth bins from the chla max depth
4. Find the max value of the pH and DO data (and respective depths) in that range
"""

import os
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from collections import OrderedDict
from shapely.geometry import Point
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
np.set_printoptions(suppress=True)


def get_season(timestamp):
    month = timestamp.month
    year = timestamp.year
    if month == 12:
        season = 'Winter'
        year = year + 1
    elif month in [1, 2]:
        season = 'Winter'
    elif month in [3, 4, 5]:
        season = 'Spring'
    elif month in [6, 7, 8]:
        season = 'Summer'
    elif month in [9, 10, 11]:
        season = 'Fall'
    else:
        season = 'Unknown'
    return f"{season} {year}"


def main(fname, project, savedir):
    # get NOAA bottom trawl survey modified strata polygons
    kwargs = dict()
    kwargs['extend'] = True
    strata_mapping = cf.return_noaa_polygons(**kwargs)

    comment = "Synthesis of chl-a profile maximum data from glider-based measurements collected during the RMI ecogliders project."
    
    # initialize dictionary to append data from glider deployments
    data = {
        "coords": {
            "time": {"dims": "time",
                     "data": np.array([], dtype='float32'),
                     "attrs": {
                         "units": "seconds since 1970-01-01T00:00:00Z",
                         "time_origin": "01-JAN-1970 00:00:00"
                     }
                }
        },
        "attrs": {
            "comment": comment,
            "acknoledgement": "This work was supported by New Jersey's Research and Monitoring Initiative (RMI) (New Jersey Department of Environmental Protection, New Jersey Board of Public Utilities).",
            "program": "An ecological and oceanographic baseline to inform offshore wind development over the continental shelf off the coast of New Jersey",
            "project": "RMI Eco-gliders"
        },
        "dims": "time",
        "data_vars": {
            "deployment": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Deployment name"
                }
            },
            "shelf_location": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Location of data point on the NES, according to modified NOAA bottom trawl survey strata. Options: inshore, midshelf, offshore"
                }
            },
            "season_year": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Year and season assignment for data point. Seasons are defined as: Winter (Dec-Feb), Spring (Mar-May), Summer (Jun-Aug), Fall (Sep-Nov)."
                }
            },
            "profile_lat": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_north",
                    "long_name": "Latitude"
                }
            },
            "profile_lon": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_east",
                    "long_name": "Longitude"
                }
            },
            "chlorophyll_a_max": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "ug L-1",
                    "long_name": "Chlorophyll-a maximum",
                    "comment": "1m depth-binned average of the data collected at the chl-a maximum of each glider profile"
                }
            },
            "chlorophyll_a_max_depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Chlorophyll-a maximum depth",
                    "comment": "Depth at which the 1m depth-binned average maximum chl-a value was recorded in the profile"
                }
            },
            "pH_near_chl_max": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "pH max near chlorophyll-a maximum",
                    "comment": "pH max within +/- 3m of the chl-a max depth"
                }
            },
            "pH_near_chl_max_depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "comment": "depth of pH max within +/- 3m of the chl-a max depth"
                }
            },
            "DO_near_chl_max": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "DO max near chlorophyll-a maximum",
                    "comment": "DO max within +/- 3m of the chl-a max depth"
                }
            },
            "DO_near_chl_max_depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "comment": "depth of DO max within +/- 3m of the chl-a max depth"
                }
            },
        },
    }

    # get glider data
    df = pd.read_csv(fname)
    
    # grab data from each deployment listed in the csv file
    for i, row in df.iterrows():
        if row['project'] == project:
            if row['use_for_analysis'] == 'yes':
                ds = xr.open_dataset(row['path_to_final_file'])
                try:
                    deployment = ds.attrs['deployment']
                except KeyError:
                    deployment = row['deployment']
                
                print(f"Processing {deployment}")

                # convert to dataframe and group by profile_time
                df = ds.to_dataframe()

                # if oxygen concentration is in umol/kg, convert to mg/L
                if 'oxygen_concentration_shifted_mgL' not in df.columns:
                    try:
                        df['oxygen_concentration_shifted_mgL'] = df['oxygen_concentration_shifted'] * 32 / 1000
                    except KeyError:
                        pass

                grouped = df.groupby('profile_time', dropna=False)
                for ts, group in grouped:
                    # check that the profile spans >10m
                    if np.nanmax(group.depth_interpolated) - np.nanmin(group.depth_interpolated) > 10:

                        # check for an MLD
                        if ~np.isnan(np.unique(group.mld)[0]):
                            # find the water column depth
                            maxdepth = np.nanmax(group.depth_interpolated)
                            glider_lon = np.nanmean(group.profile_lon)
                            glider_lat = np.nanmean(group.profile_lat)
                        
                            pdf = cf.depth_bin(group, depth_var='depth_interpolated', depth_min=0, depth_max=maxdepth, stride=1)

                            if len(pdf) > 0:
                                # determine the shelf location
                                glider_point = Point(glider_lon, glider_lat)
                                shelf_location = 'unknown'

                                for region, info in strata_mapping.items():
                                    if info['poly'].contains(glider_point):
                                        shelf_location = region
                                        break

                                data['coords']['time']['data'] = np.append(data['coords']['time']['data'], pd.to_datetime(ts).timestamp())

                                season_year = get_season(ts)

                                data['data_vars']['deployment']['data'] = np.append(data['data_vars']['deployment']['data'], deployment)
                                data['data_vars']['shelf_location']['data'] = np.append(data['data_vars']['shelf_location']['data'], shelf_location)
                                data['data_vars']['season_year']['data'] = np.append(data['data_vars']['season_year']['data'], season_year)
                                for v in ['profile_lat', 'profile_lon']:
                                    avg = np.nanmean(pdf[v].values)
                                    data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], avg)
                                        
                                pdf.rename_axis('depth_bin', inplace=True)
                                pdf.reset_index(inplace=True)
                                chlmax = np.nanmax(pdf['chlorophyll_a'])
                                data['data_vars']['chlorophyll_a_max']['data'] = np.append(data['data_vars']['chlorophyll_a_max']['data'], chlmax)
                                if np.isnan(chlmax):
                                    chlmax_depth = np.nan
                                    phmax = np.nan
                                    phmax_depth = np.nan
                                    domax = np.nan
                                    domax_depth = np.nan
                                else:
                                    chlmax_index = pdf.loc[pdf['chlorophyll_a'] == chlmax].index[0]
                                    chlmax_depth = pdf.iloc[chlmax_index].depth_interpolated

                                    # subset the dataframe to +/- 3m depth bins of the chl-a max depth
                                    chlmax_df = pdf.iloc[max(0, chlmax_index - 3):chlmax_index + 4]
                                    
                                    # get pH
                                    try:
                                        phmax = np.nanmax(chlmax_df['pH'])
                                    except KeyError:
                                        phmax = np.nan
                                    if np.isnan(phmax):
                                        phmax_depth = np.nan
                                    else:
                                        phmax_depth = chlmax_df.loc[chlmax_df['pH'] == phmax].depth_interpolated
                                    
                                    # get DO
                                    try:
                                        domax = np.nanmax(chlmax_df['oxygen_concentration_shifted_mgL'])
                                    except KeyError:
                                        domax = np.nan
                                    if np.isnan(domax):
                                        domax_depth = np.nan
                                    else:
                                        domax_depth = chlmax_df.loc[chlmax_df['oxygen_concentration_shifted_mgL'] == domax].depth_interpolated
                                    
                                # add variables to the dataset
                                data['data_vars']['chlorophyll_a_max_depth']['data'] = np.append(data['data_vars']['chlorophyll_a_max_depth']['data'], chlmax_depth)
                                data['data_vars']['pH_near_chl_max']['data'] = np.append(data['data_vars']['pH_near_chl_max']['data'], phmax)
                                data['data_vars']['pH_near_chl_max_depth']['data'] = np.append(data['data_vars']['pH_near_chl_max_depth']['data'], phmax_depth)
                                data['data_vars']['DO_near_chl_max']['data'] = np.append(data['data_vars']['DO_near_chl_max']['data'], domax)
                                data['data_vars']['DO_near_chl_max_depth']['data'] = np.append(data['data_vars']['DO_near_chl_max_depth']['data'], domax_depth)

                print(f'Finished {deployment}')

    # save data as netcdf
    ds = xr.Dataset.from_dict(data)

    # add created time to global attrs
    datetime_format = '%Y-%m-%dT%H:%M:%SZ'
    created = dt.datetime.now(dt.UTC).strftime(datetime_format)  # creation time Timestamp
    time_start = dt.datetime.fromtimestamp(np.nanmin(ds.time.values), dt.UTC).strftime('%Y-%m-%d')
    time_end = dt.datetime.fromtimestamp(np.nanmax(ds.time.values), dt.UTC).strftime('%Y-%m-%d')
    start_yr = dt.datetime.fromtimestamp(np.nanmin(ds.time.values), dt.UTC).strftime('%Y')
    end_yr = dt.datetime.fromtimestamp(np.nanmax(ds.time.values), dt.UTC).strftime('%Y')

    global_attributes = OrderedDict([
        ('date_created', created),
        ('date_modified', created),
        ('time_coverage_start', time_start),
        ('time_coverage_end', time_end),
        ('creator_email', 'lgarzio@marine.rutgers.edu'),
        ('creator_name', 'Lori Garzio'),
        ('creator_url', 'rucool.marine.rutgers.edu'),
        ('institution', 'Rutgers University'),
        ('contributor_name', 'Grace Saba,Lori Garzio'),
        ('contributor_role', 'Principal Investigator,Data Management')
    ])

    global_attributes.update(ds.attrs)

    ds = ds.assign_attrs(global_attributes)
    ds = ds.sortby(ds.time)

    # Add compression to all variables
    encoding = {}
    for k in ds.data_vars:
        if k not in ['deployment', 'shelf_location', 'season_year']:
            encoding[k] = {'zlib': True, 'complevel': 1}

    encoding['time'] = dict(zlib=False, _FillValue=False, dtype=np.double)

    save_file = os.path.join(savedir, f'RMI_gliderdata_chlamax_{start_yr}_{end_yr}.nc')
    ds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    csv_summary = '/Users/garzio/Documents/rucool/Saba/RMI/glider_deployments.csv'
    project = 'RMI'
    savedir = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/data'
    main(csv_summary, project, savedir)
