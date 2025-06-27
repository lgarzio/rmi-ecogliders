#!/usr/bin/env python

"""
Author: Lori Garzio on 8/10/2021
Last modified: 6/24/2025
Quickly plot xsections of glider data variables, adding mixed layer depth if specified
"""

import os
import xarray as xr
import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
import cmocean as cmo
import functions.common as cf
import functions.plotting as pf
plt.rcParams.update({'font.size': 13})


def main(fname, add_mld):
    savedir = os.path.join(os.path.dirname(fname), 'xsections')
    if add_mld:
        savedir = os.path.join(os.path.dirname(fname), 'xsections_mld')
    os.makedirs(savedir, exist_ok=True)

    ds = xr.open_dataset(fname)
    try:
        deploy = ds.attrs['deployment']
    except KeyError:
        f = fname.split('/')[-1]
        deploy = f'{f.split("-")[0]}-{f.split("-")[1]}'  # get the deployment from the filename

    t0str = pd.to_datetime(np.nanmin(ds.time)).strftime('%Y-%m-%dT%H:%M')
    t1str = pd.to_datetime(np.nanmax(ds.time)).strftime('%Y-%m-%dT%H:%M')

    root_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(root_dir)
    configdir = os.path.join(parent_dir, 'config')

    with open(os.path.join(configdir, 'plot_vars.yml')) as f:
        plt_vars = yaml.safe_load(f)

    for pv, info in plt_vars.items():
        try:
            variable = ds[pv]
        except KeyError:
            continue

        # plot xsection
        if len(variable) > 1:
            fig, ax = plt.subplots(figsize=(12, 6))
            plt.subplots_adjust(left=0.1)
            figttl_xsection = f'{deploy} {variable.attrs['long_name']}\n{t0str} to {t1str}'
            clab = f'{variable.attrs['long_name']} ({variable.attrs['units']})'
            xargs = dict()
            xargs['clabel'] = clab
            xargs['title'] = figttl_xsection
            xargs['date_fmt'] = '%m-%d'
            xargs['grid'] = True
            xargs['cmap'] = info['cmap']
            pf.xsection(fig, ax, ds.time.values, ds.depth_interpolated.values, variable.values, **xargs)

            if add_mld:
                ax.plot(ds.time.values, ds.mld.values, color='black', linewidth=1.5, ls='-')
            
            sfilename = f'{deploy}_xsection_{pv}.png'
            sfile = os.path.join(savedir, sfilename)
            plt.savefig(sfile, dpi=300)
            plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2023/ru39-20231103T1413/ncei_pH/ru39-20231103T1413-delayed_mld.nc'
    add_mld = True  # add MLD to the plot, True or False
    main(ncfile, add_mld)
