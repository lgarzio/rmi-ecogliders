#!/usr/bin/env python

"""
Author: Lori Garzio on 7/1/2025
Last modified: 7/8/2025
Glider data wrangler, export as NetCDF. The netcdf files have already been QC'd and processed, so this script
is for synthesizing the surface, bottom, and maximum (chl, DO, pH) data from glider deployments and grouping them 
according to region (inshore, midshelf, offshore) based on NOAA strata.
Surface data: for each profile, grabbed data between 1-5m depth took the average of that and appended to the dataset.
Bottom data: for each profile, found the max glider depth, compared that to GEBCO bathymetry file using profile lat/lon. 
If the max glider depth was within +/- 20% of the depth from the bathymetry file, grabbed the bottom 4m of data from that
profile, took the average of that and appended to the dataset.
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


def main(fname, project, surf_bot, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-75.5, -71.5, 38, 41]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))
    
    # get NOAA bottom trawl survey modified strata polygons
    kwargs = dict()
    kwargs['extend'] = True
    strata_mapping = cf.return_noaa_polygons(**kwargs)

    if surf_bot == 'surface':
        comment = "Synthesis of surface data (average of profile data from 1-5m depth) from glider-based measurements collected during the RMI ecogliders project."
    elif surf_bot == 'bottom':
        comment = "Synthesis of bottom data (average of data from the bottom 4m of each profile) from glider-based measurements collected during the RMI ecogliders project."
    elif surf_bot == 'maximum':
        comment = "Synthesis of profile maximum data from glider-based measurements collected during the RMI ecogliders project."
    else:
        raise ValueError("surf_bot must be either surface, bottom or maximum")
    
    # initialize dictionary to append surface or bottom data from glider deployments
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
            "depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Depth",
                    "description": "Average depth of glider sample collection"
                }
            },
            "bottom_depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Bottom Depth",
                    "description": "Depth of water column at the sampling location",
                    "comment": "Bottom depth was determined using the profile coordinates and GEBCO’s gridded global "
                               "bathymetry (https://download.gebco.net/)"
                }
            },
        },
    }

    # add additional science variables to dataset
    root_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(root_dir)
    configdir = os.path.join(parent_dir, 'config')

    with open(os.path.join(configdir, 'sci_vars.yml')) as f:
        sci_vars = yaml.safe_load(f)
    
    if surf_bot == 'maximum':
        # remove "depth" from the data_vars dict
        del data['data_vars']['depth']
        
        # add maximum data variables and the depth of the maxima to dataset
        for sv in ['chlorophyll_a', 'oxygen_concentration_shifted_mgL', 'pH']:
            data['data_vars'][f'{sv}_max'] = {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": sci_vars[sv]
            }
            data['data_vars'][f'{sv}_max']['attrs']['comment'] = f"1m depth-binned average of the data collected at the {sv} {surf_bot} of each glider profile"

            attrs = dict(
                units="m",
                long_name=f"{sci_vars[sv]['long_name']} maximum depth",
                comment=f"Depth at which the 1m depth-binned average maximum {sv} value was recorded in the profile",
            )
            data['data_vars'][f'{sv}_max_depth'] = {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": attrs
            }
        
    else:
        for sv, vals in sci_vars.items():
            data['data_vars'][sv] = {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": vals
            }
            data['data_vars'][sv]['attrs']['comment'] = f"Average of the data collected at the {surf_bot} of each glider profile"

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
                
                grouped = df.groupby('profile_time', dropna=False)
                for ts, group in grouped:
                    # check that the profile spans >10m
                    if np.nanmax(group.depth_interpolated) - np.nanmin(group.depth_interpolated) > 10:
                        # find the water column depth
                        maxdepth = np.nanmax(group.depth_interpolated)
                        glider_lon = np.nanmean(group.profile_lon)
                        glider_lat = np.nanmean(group.profile_lat)
                        lat_idx = abs(bathy.lat.values - glider_lat).argmin()
                        lon_idx = abs(bathy.lon.values - glider_lon).argmin()
                        station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                        if surf_bot == 'surface':
                            # grab top 1-5m of data
                            pdf = group.loc[np.logical_and(group.depth_interpolated <= 5, group.depth_interpolated >= 1)]
                        elif surf_bot == 'bottom':
                            # compare the glider depth to the global bathymetry file
                            # if the glider sample is within +/- 20% of the water column,
                            # grab the bottom 4m of data
                            depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                            if np.logical_and(maxdepth > depth_threshold[0], maxdepth < depth_threshold[1]):
                                pdf = group.loc[group.depth_interpolated >= (maxdepth - 4)]
                            else:
                                pdf = pd.DataFrame()
                        elif surf_bot == 'maximum':
                            # bin the data into 1m depth bins averages, and grab the maximum of the profile
                            pdf = cf.depth_bin(group, depth_var='depth_interpolated', depth_min=0, depth_max=maxdepth, stride=1)
                        else:
                            raise ValueError("surf_bot must be either surface, bottom or maximum")

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

                            for v in data['data_vars']:
                                if v == 'bottom_depth':
                                    # add bottom depth to the dataset
                                    data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'],
                                                                            station_water_depth)
                                elif v == 'deployment':
                                    # add deployment name to the dataset
                                    data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'],
                                                                            deployment)
                                elif v == 'shelf_location':
                                    # add shelf location to the dataset
                                    data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], shelf_location)
                                elif v == 'season_year':
                                    data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], season_year)
                                else:
                                    try:
                                        avg = np.nanmean(pdf[v].values)
                                        data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], avg)
                                    except KeyError:
                                        if surf_bot == 'maximum':
                                            if v.endswith('_max'):
                                                # grab the maximum value from the profile
                                                tempv = v.replace('_max', '')
                                                try:
                                                    max_val = np.nanmax(pdf[tempv])
                                                except KeyError:
                                                    if tempv == 'oxygen_concentration_shifted_mgL':
                                                        try:
                                                            max_val = np.nanmax(pdf['oxygen_concentration_shifted'] * 32 / 1000)
                                                        except KeyError:
                                                            max_val = np.nan
                                                    else:
                                                        max_val = np.nan
                                                data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], max_val)
                                            elif v.endswith('_max_depth'):
                                                # grab the depth of the maximum value
                                                tempv = v.replace('_max_depth', '')

                                                # figure out if the variable is there
                                                try:
                                                    tempv_data = pdf[tempv]
                                                except KeyError:
                                                    if tempv == 'oxygen_concentration_shifted_mgL':
                                                        try:
                                                            tempv_data = pdf['oxygen_concentration_shifted']
                                                            tempv = 'oxygen_concentration_shifted'
                                                        except KeyError:
                                                            tempv_data = np.nan
                                                    else:
                                                        tempv_data = np.nan
                                                
                                                # if the maximum value is nan, then the depth is also nan
                                                if np.isnan(np.nanmax(tempv_data)):
                                                    vmax_depth = np.nan
                                               
                                               # otherwise, find the depth at the variable maximum
                                                else:
                                                    try:
                                                        vmax_depth = pdf.loc[pdf[tempv] == np.nanmax(pdf[tempv])].depth_interpolated.values[0]
                                                    except KeyError:
                                                        vmax_depth = np.nan
                                                data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], vmax_depth)
                                        else:
                                            if v == 'oxygen_concentration_shifted_mgL':
                                                try:
                                                    avg = np.nanmean(pdf['oxygen_concentration_shifted'].values * 32 / 1000)  # convert from umol/L to mg/L
                                                    data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], avg)
                                                except KeyError:
                                                    data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], np.nan)
                                            else:
                                                data['data_vars'][v]['data'] = np.append(data['data_vars'][v]['data'], np.nan)
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

    save_file = os.path.join(savedir, f'RMI_gliderdata_{surf_bot}_{start_yr}_{end_yr}.nc')
    ds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    csv_summary = '/Users/garzio/Documents/rucool/Saba/RMI/glider_deployments.csv'
    project = 'RMI'
    surf_bot = 'maximum'  # grab surface, bottom or maximum data
    savedir = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/data'
    main(csv_summary, project, surf_bot, savedir)
