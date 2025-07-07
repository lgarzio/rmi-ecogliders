#!/usr/bin/env python

"""
Author: Lori Garzio on 6/27/2025
Last modified: 6/27/2025
Create 2-panel profile plots for each profile for multiple variables
"""

import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import functions.common as cf
plt.rcParams.update({'font.size': 12})
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


def main(fname, mldvers):
    savedir = os.path.join(os.path.dirname(fname), f'profiles_mld_{mldvers}')
    os.makedirs(savedir, exist_ok=True)
    
    ds = xr.open_dataset(fname)
    ds = ds.sortby(ds.time)
    try:
        deploy = ds.attrs['deployment']
    except KeyError:
        f = fname.split('/')[-1]
        deploy = f'{f.split("-")[0]}-{f.split("-")[1]}'  # get the deployment from the filename

    df = ds.to_dataframe()
    grouped = df.groupby('profile_time')

    for group in grouped:
        if np.sum(~np.isnan(group[1]['mld'])) > 0:
            tstr = group[0].strftime("%Y-%m-%dT%H%M%SZ")
            temp = group[1]['temperature']
            sal = group[1]['salinity']
            dens = group[1]['density']
            chla = group[1]['chlorophyll_a']
            pH = group[1]['pH']

            fig, axs = plt.subplots(1, 2, figsize=(15, 11), sharey=True)
            ms = 20
            
            # Main axis for temperature
            axs[0].scatter(temp, group[1]['depth_interpolated'], s=ms, color='red')
            axs[0].set_xlabel('Temperature (°C)', color='red')
            axs[0].tick_params(axis='x', colors='red')
            axs[0].set_ylabel('Depth (m)')

            # Add second x-axis for salinity
            ax2 = axs[0].twiny()
            ax2.scatter(sal, group[1]['depth_interpolated'], s=ms, color='darkturquoise')
            ax2.set_xlabel('Salinity', color='darkturquoise')
            ax2.tick_params(axis='x', colors='darkturquoise')

            # Add third x-axis for density
            ax3 = axs[0].twiny()
            ax3.scatter(dens, group[1]['depth_interpolated'], s=ms, color='indigo')
            ax3.spines['top'].set_position(('outward', 40))  # Offset the axis
            ax3.set_xlabel('Density', color='indigo')
            ax3.tick_params(axis='x', colors='indigo')

            # Add chlorophyll-a on second plot
            axs[1].scatter(chla, group[1]['depth_interpolated'], s=ms, color='green')
            axs[1].set_xlabel('Chlorophyll-a', color='green')
            axs[1].tick_params(axis='x', colors='green')

            # Add pH on the plot with chla
            ax5 = axs[1].twiny()
            ax5.scatter(pH, group[1]['depth_interpolated'], s=ms, color='darkorange')
            ax5.set_xlabel('pH', color='darkorange')
            ax5.tick_params(axis='x', colors='darkorange')

            axs[0].invert_yaxis()
            axs[0].set_ylabel('Depth (m)')

            fig.suptitle(f'{deploy.split("-")[0]} {tstr}')

            if mldvers == 'n2':
                axs[0].axhline(y=np.unique(group[1]['mld']), ls='--', c='k')
                axs[1].axhline(y=np.unique(group[1]['mld']), ls='--', c='k', label='MLD')
                axs[1].legend(fontsize=10)
            elif mldvers == '2way':
                axs[0].axhline(y=np.unique(group[1]['mld_upper']), ls='--', c='k')
                axs[0].axhline(y=np.unique(group[1]['mld_lower']), ls='--', c='k')
                axs[1].axhline(y=np.unique(group[1]['mld_upper']), ls='--', c='k')
                axs[1].axhline(y=np.unique(group[1]['mld_lower']), ls='--', c='k')

            plt.subplots_adjust(top=.86, bottom=0.06)

            sfile = os.path.join(savedir, f'{tstr}_profile.png')
            plt.savefig(sfile, dpi=300)
            plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2023/ru39-20230817T1520/ncei_pH/ru39-20230817T1520-delayed_mld2way.nc'
    mld_version = '2way'  # 'n2' '2way
    main(ncfile, mld_version)
