#!/usr/bin/env python

"""
Author: Lori Garzio on 6/25/2025
Last modified: 6/25/2025
Plot one glider track
"""

import os
import xarray as xr
import pandas as pd
import cool_maps.plot as cplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import functions.common as cf
plt.rcParams.update({'font.size': 13})


def main(fname, savedir, add_strata):
    ds = xr.open_dataset(fname)

    try:
        deploy = ds.attrs['deployment']
    except KeyError:
        f = fname.split('/')[-1]
        deploy = f'{f.split("-")[0]}-{f.split("-")[1]}'  # get the deployment from the filename

    df = pd.DataFrame({'lon': ds.profile_lon.values, 'lat': ds.profile_lat.values})
    df = df.drop_duplicates()
    
    extent = [-75.5, -71.5, 38, 41]
    
    kwargs = dict()
    kwargs['figsize'] = (9, 8)
    kwargs['coast'] = 'low'  # low full
    kwargs['oceancolor'] = 'none'
    kwargs['decimal_degrees'] = True
    kwargs['bathymetry'] = True
    #kwargs['bathymetry_file'] = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    kwargs['bathymetry_method'] = 'topo_log'
    fig, ax = cplt.create(extent, **kwargs)

    ax.scatter(df['lon'], df['lat'], color='magenta', marker='.', s=6, transform=ccrs.PlateCarree(), zorder=30)
    
    sfilename = f'{deploy}_glider_track.png'

    if add_strata:
        kwargs = dict()
        kwargs['extend'] = True
        strata_mapping = cf.return_noaa_polygons(**kwargs)
        for key, values in strata_mapping.items():
            outside_poly = values['poly']
            xx, yy = outside_poly.exterior.xy
            ax.plot(xx, yy, color=values['color'], lw=4, transform=ccrs.PlateCarree(), zorder=20)
        sfilename = f'{deploy}_glider_track_strata.png'
    
    sfile = os.path.join(savedir, sfilename)

    plt.subplots_adjust(top=.92, bottom=0.08, right=1, left=0)
    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    fname = '/Users/garzio/Documents/rucool/Saba/gliderdata/2023/ru40-20231103T1421/ncei_dmon/ru40-20231103T1421-20231115T1612-delayed_mld_combined.nc'
    savedir = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/maps'
    add_strata = True  # add NOAA strata to the plot, True or False
    main(fname, savedir, add_strata)
