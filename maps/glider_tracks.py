#!/usr/bin/env python

"""
Author: Lori Garzio on 6/25/2025
Last modified: 6/25/2025
Plot all RMI glider tracks
"""

import os
import xarray as xr
import pandas as pd
import cool_maps.plot as cplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import functions.common as cf
plt.rcParams.update({'font.size': 13})


def main(fname, project, savedir, add_strata):
    df = pd.read_csv(fname)

    data = dict()

    # grab data from each deployment listed in the csv file
    for i, row in df.iterrows():
        if row['project'] == project:
            if row['use_for_analysis'] == 'yes':
                ds = xr.open_dataset(row['path_to_final_file'])
                try:
                    deploy = ds.attrs['deployment']
                except KeyError:
                    deploy = row['deployment']

                coordf = pd.DataFrame({'lon': ds.profile_lon.values, 'lat': ds.profile_lat.values})
                coordf = coordf.drop_duplicates()

                data[deploy] = dict()
                data[deploy]['lon'] = coordf.lon.values
                data[deploy]['lat'] = coordf.lat.values

    extent = [-75.5, -71.5, 38, 41]
    
    kwargs = dict()
    kwargs['figsize'] = (9, 8)
    kwargs['coast'] = 'full'  # low full
    kwargs['oceancolor'] = 'none'
    kwargs['decimal_degrees'] = True
    kwargs['bathymetry'] = True
    #kwargs['bathymetry_file'] = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    kwargs['bathymetry_method'] = 'topo_log'
    fig, ax = cplt.create(extent, **kwargs)

    for deploy, coords in data.items():
        ax.scatter(coords['lon'], coords['lat'], color='magenta', marker='.', s=6, 
                   transform=ccrs.PlateCarree(), zorder=30)
    
    sfilename = 'RMI_glider_tracks.png'

    if add_strata:
        kwargs = dict()
        kwargs['extend'] = True
        strata_mapping = cf.return_noaa_polygons(**kwargs)
        for key, values in strata_mapping.items():
            outside_poly = values['poly']
            xx, yy = outside_poly.exterior.xy
            ax.plot(xx, yy, color=values['color'], lw=4, transform=ccrs.PlateCarree(), zorder=20)
        sfilename = 'RMI_glider_tracks_strata.png'
    
    sfile = os.path.join(savedir, sfilename)

    plt.subplots_adjust(top=.92, bottom=0.08, right=1, left=0)
    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    csv_summary = '/Users/garzio/Documents/rucool/Saba/RMI/glider_deployments.csv'
    project = 'RMI'
    savedir = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/maps'
    add_strata = True  # add NOAA strata to the plot, True or False
    main(csv_summary, project, savedir, add_strata)
