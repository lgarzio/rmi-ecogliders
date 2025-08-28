#!/usr/bin/env python

"""
Author: Lori Garzio on 6/25/2025
Last modified: 8/28/2025
Plot all RMI glider tracks
"""

import os
import xarray as xr
import pandas as pd
import cool_maps.plot as cplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import functions.common as cf
import functions.plotting as pf
plt.rcParams.update({'font.size': 13})


def main(fname, project, savedir, add_strata, add_wea, extent):
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
    
    kwargs = dict()
    kwargs['figsize'] = (9, 8)
    kwargs['coast'] = 'full'  # low full
    kwargs['oceancolor'] = 'none'
    kwargs['decimal_degrees'] = True
    kwargs['bathymetry'] = True
    #kwargs['bathymetry_file'] = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    kwargs['bathymetry_method'] = 'topo_log'
    fig, ax = cplt.create(extent, **kwargs)

    if add_wea:
        lease = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/boem-renewable-energy-shapefiles_0/Offshore_Wind_Leases_outlines.shp'
        kwargs = dict()
        kwargs['edgecolor'] = 'gray'
        kwargs['facecolor'] = '#afb0ae'  # light gray
        kwargs['zorder'] = 10
        kwargs['alpha'] = 0.8
        pf.map_add_shapefile_outlines(ax, lease, **kwargs)

    for deploy, coords in data.items():
        ax.scatter(coords['lon'], coords['lat'], color='magenta', marker='.', s=6, 
                   transform=ccrs.PlateCarree(), zorder=30)
    
    sfilename = f'{project}_glider_tracks.png'

    if add_strata:
        kwargs = dict()
        kwargs['extend'] = True
        strata_mapping = cf.return_noaa_polygons(**kwargs)
        for key, values in strata_mapping.items():
            outside_poly = values['poly']
            xx, yy = outside_poly.exterior.xy
            ax.plot(xx, yy, color=values['color'], lw=4, transform=ccrs.PlateCarree(), zorder=values['zorder'],
                    label=f'{key} stratum')
        sfilename = f'{project}_glider_tracks_strata.png'
    
    sfile = os.path.join(savedir, sfilename)

    plt.subplots_adjust(top=.92, bottom=0.08, right=1, left=0)
    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    csv_summary = '/Users/garzio/Documents/rucool/Saba/RMI/glider_deployments.csv'
    project = 'RMI'
    savedir = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/maps'
    add_strata = True  # add NOAA strata to the plot, True or False
    add_wea = True  # add Wind Energy Areas to the plot, True or False
    extent = [-75, -72, 38.4, 40.7]  # [-75.5, -71.5, 38, 41]  [-75, -72, 38.4, 40.8]
    main(csv_summary, project, savedir, add_strata, add_wea, extent)
