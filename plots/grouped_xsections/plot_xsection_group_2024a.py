#!/usr/bin/env python

"""
Author: Lori Garzio on 6/27/2025
Last modified: 7/15/2025
Plot cross-sections of data from groups of glider deployments with shared axes.

"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean as cmo
import functions.plotting as pf
import functions.oxy_colormap_mods as ocm
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console
plt.rcParams.update({'font.size': 17})
#plt.rcParams.update({'font.size': 15})


def main(savedir):
    os.makedirs(savedir, exist_ok=True)

    spring_ph = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru39-20240429T1522/ncei_pH/ru39-20240429T1522-delayed_mld.nc')
    spring_do = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru40-20240429T1528/ncei_dmon/ru40-20240429T1528-delayed_mld.nc')
    latespring =  xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru43-20240612T1658/ncei_pH/ru43-20240612T1658-delayed_mld.nc') 
    summer_ph = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru39-20240723T1442/ncei_pH/ru39-20240723T1442-delayed_mld.nc')
    summer_do = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru40-20240723T1600/ncei_dmon/ru40-20240723T1600-delayed_mld.nc')
    #latesummer = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru43-20240904T1539/ncei_pH/ru43-20240904T1539-delayed_mld.nc')
    #fall_ph = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru43-20240904T1539/ncei_pH/ru43-20240904T1539-delayed_mld.nc')
    #fall_do = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2024/ru40-20241021T1654/ncei_dmon/ru40-20241021T1654-delayed_mld.nc')
    #winter_ph = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2025/ru39-20250226T1700/ncei_pH/ru39-20250226T1700-delayed_mld.nc')
    #winter_do = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2025/ru40-20250226T1704/ncei_dmon/ru40-20250226T1704-delayed_mld.nc')

    # # make summary dataframe to export as csv
    # summary_dict = dict(variable=['temperature', 'chl', 'DO_mgL', 'pH', 'omega'],
    #                     gapfill_min=[],
    #                     gapfill_max=[],
    #                     fall_min=[],
    #                     fall_max=[])

    # plot cross sections of all deployments: temperature, chl, DO, pH, omega
    fig, axs = plt.subplots(6, 3, figsize=(16, 18), sharex='col', sharey=True)
    plt.subplots_adjust(top=.96, bottom=0.08, right=.95, left=.06, hspace=0.05, wspace=0.05)
    ax1 = axs[0, 0]  # spring temperature
    ax2 = axs[0, 1]  # latespring temperature
    ax3 = axs[0, 2]  # summer temperature
    ax5 = axs[1, 0]  # spring salinity
    ax6 = axs[1, 1]  # latespring salinity
    ax7 = axs[1, 2]  # summer salinity
    ax9 = axs[2, 0]  # spring chl
    ax10 = axs[2, 1]  # latespring chl
    ax11 = axs[2, 2]  # summer chl
    ax13 = axs[3, 0]  # spring DO - not deployed
    ax14 = axs[3, 1]  # latespring DO
    ax15 = axs[3, 2]  # summer DO
    ax17 = axs[4, 0]  #  spring pH
    ax18 = axs[4, 1]  #  latespring pH
    ax19 = axs[4, 2]  #  summer pH
    ax20 = axs[5, 0]  # spring omega
    ax21 = axs[5, 1]   # latespring omega
    ax22 = axs[5, 2]  # summer omega

    # plot temperature
    kwargs = dict()
    kwargs['clabel'] = 'Temperature ({})'.format(r'$\rm ^oC$')
    kwargs['cmap'] = cmo.cm.thermal
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [7, 23]
    kwargs['title'] = 'Spring 2024'
    kwargs['ylims'] = [-4, 105]

    pf.xsection(fig, ax1, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.temperature.values, **kwargs)  # spring temperature

    kwargs['title'] = 'Early Summer 2024'
    kwargs['ylabel'] = None
    pf.xsection(fig, ax2, latespring.time.values, latespring.depth_interpolated.values, latespring.temperature.values, **kwargs)  # latespring temperature

    kwargs['title'] = 'Summer 2024'
    kwargs['ylabel'] = None
    kwargs['colorbar'] = True
    pf.xsection(fig, ax3, summer_ph.time.values, summer_ph.depth_interpolated.values, summer_ph.temperature.values, **kwargs)  # summer temperature

    # print('temperature min/max')
    # print(f'gapfill {np.nanmin(ds_gapfill.temperature.values)} / {np.nanmax(ds_gapfill.temperature.values)}')
    # print(f'fall {np.nanmin(ds_fall_ph.temperature.values)} / {np.nanmax(ds_fall_ph.temperature.values)}')
    # summary_dict['gapfill_min'].append(np.round(np.nanmin(ds_gapfill.temperature.values), 2))
    # summary_dict['gapfill_max'].append(np.round(np.nanmax(ds_gapfill.temperature.values), 2))
    # summary_dict['fall_min'].append(np.round(np.nanmin(ds_fall_ph.temperature.values), 2))
    # summary_dict['fall_max'].append(np.round(np.nanmax(ds_fall_ph.temperature.values), 2))

    # plot salinity
    kwargs = dict()
    kwargs['clabel'] = 'Salinity (PSU)'
    kwargs['cmap'] = cmo.cm.haline
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [29, 35]

    pf.xsection(fig, ax5, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.salinity.values, **kwargs)  # spring chl

    kwargs['ylabel'] = None
    pf.xsection(fig, ax6, latespring.time.values, latespring.depth_interpolated.values, latespring.salinity.values, **kwargs)  # latespring chl
    
    kwargs['colorbar'] = True
    pf.xsection(fig, ax7, summer_ph.time.values, summer_ph.depth_interpolated.values, summer_ph.salinity.values, **kwargs)  # summer chl

    # plot chl
    kwargs = dict()
    kwargs['clabel'] = 'Chlorophyll a (ug/L)'
    kwargs['cmap'] = cmo.cm.algae
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [2, 10]

    pf.xsection(fig, ax9, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.chlorophyll_a.values, **kwargs)  # spring chl

    kwargs['ylabel'] = None
    pf.xsection(fig, ax10, latespring.time.values, latespring.depth_interpolated.values, latespring.chlorophyll_a.values, **kwargs)  # latespring chl
    
    kwargs['colorbar'] = True
    pf.xsection(fig, ax11, summer_ph.time.values, summer_ph.depth_interpolated.values, summer_ph.chlorophyll_a.values, **kwargs)  # summer chl

    # print('chl min/max')
    # print(f'gapfill {np.nanmin(ds_gapfill.chlorophyll_a.values)} / {np.nanmax(ds_gapfill.chlorophyll_a.values)}')
    # print(f'fall {np.nanmin(ds_fall_ph.chlorophyll_a.values)} / {np.nanmax(ds_fall_ph.chlorophyll_a.values)}')

    # summary_dict['gapfill_min'].append(np.round(np.nanmin(ds_gapfill.chlorophyll_a.values), 2))
    # summary_dict['gapfill_max'].append(np.round(np.nanmax(ds_gapfill.chlorophyll_a.values), 2))
    # summary_dict['fall_min'].append(np.round(np.nanmin(ds_fall_ph.chlorophyll_a.values), 2))
    # summary_dict['fall_max'].append(np.round(np.nanmax(ds_fall_ph.chlorophyll_a.values), 2))

    # DO
    # get color map
    kwargs = dict()
    kwargs['breaks'] = [3, 5]
    kwargs['blue'] = False
    mymap = ocm.cm_partialturbo_r(**kwargs)
    #mymap = ocm.cm_rogg()

    kwargs = dict()
    kwargs['clabel'] = 'Oxygen (mg/L)'
    kwargs['title'] = None
    kwargs['cmap'] = mymap
    kwargs['xlabel'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [2, 9]
    kwargs['date_fmt'] = '%m-%d'

    pf.xsection(fig, ax13, spring_do.time.values, spring_do.depth_interpolated.values, spring_do.oxygen_concentration_shifted_mgL.values, **kwargs)  # spring DO

    kwargs['ylabel'] = None
    kwargs['xticklabels'] = None
    #lsdo_mgL = latespring.oxygen_concentration_shifted.values * 32 / 1000  # convert from umol/L to mg/L
    pf.xsection(fig, ax14, latespring.time.values, latespring.depth_interpolated.values, latespring.oxygen_concentration_shifted_mgL.values, **kwargs)  # latespring DO

    kwargs['colorbar'] = True
    #summerdo_mgL = summer_do.oxygen_concentration_shifted.values * 32 / 1000  # convert from umol/L to mg/L
    pf.xsection(fig, ax15, summer_do.time.values, summer_do.depth_interpolated.values, summer_do.oxygen_concentration_shifted_mgL.values, **kwargs)  # summer DO

    # print('DO min/max')
    # print(f'gapfill: none')
    # print(f'fall {np.nanmin(ds_fall_do.oxygen_concentration_shifted_mgL.values)} / {np.nanmax(ds_fall_do.oxygen_concentration_shifted_mgL.values)}')

    # summary_dict['gapfill_min'].append('nan')
    # summary_dict['gapfill_max'].append('nan')
    # summary_dict['fall_min'].append(np.round(np.nanmin(ds_fall_do.oxygen_concentration_shifted_mgL.values), 2))
    # summary_dict['fall_max'].append(np.round(np.nanmax(ds_fall_do.oxygen_concentration_shifted_mgL.values), 2))

    # pH
    kwargs = dict()
    kwargs['clabel'] = 'pH'
    kwargs['title'] = None
    kwargs['cmap'] = cmo.cm.matter
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [7.6, 8.1]

    pf.xsection(fig, ax17, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.pH.values, **kwargs)  # spring pH

    kwargs['ylabel'] = None
    pf.xsection(fig, ax18, latespring.time.values, latespring.depth_interpolated.values, latespring.pH.values, **kwargs)  # latespring pH

    kwargs['colorbar'] = True
    pf.xsection(fig, ax19, summer_ph.time.values, summer_ph.depth_interpolated.values, summer_ph.pH.values, **kwargs)  # summer pH

    # print('pH min/max')
    # print(f'gapfill {np.nanmin(ds_gapfill.pH.values)} / {np.nanmax(ds_gapfill.pH.values)}')
    # print(f'fall {np.nanmin(ds_fall_ph.pH.values)} / {np.nanmax(ds_fall_ph.pH.values)}')

    # summary_dict['gapfill_min'].append(np.round(np.nanmin(ds_gapfill.pH.values), 2))
    # summary_dict['gapfill_max'].append(np.round(np.nanmax(ds_gapfill.pH.values), 2))
    # summary_dict['fall_min'].append(np.round(np.nanmin(ds_fall_ph.pH.values), 2))
    # summary_dict['fall_max'].append(np.round(np.nanmax(ds_fall_ph.pH.values), 2))

    # omega
    kwargs = dict()
    kwargs['clabel'] = 'Omega'
    kwargs['title'] = None
    kwargs['cmap'] = cmo.cm.matter
    kwargs['xlabel'] = None
    kwargs['colorbar'] = None
    #kwargs['vlims'] = [1, 3]
    kwargs['vlims'] = [0.7, 3]
    kwargs['date_fmt'] = '%m-%d'

    pf.xsection(fig, ax20, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.aragonite_saturation_state.values, **kwargs)  # spring omega

    kwargs['ylabel'] = None
    pf.xsection(fig, ax21, latespring.time.values, latespring.depth_interpolated.values, latespring.aragonite_saturation_state.values, **kwargs)  # latespring omega

    kwargs['colorbar'] = True
    pf.xsection(fig, ax22, summer_ph.time.values, summer_ph.depth_interpolated.values, summer_ph.aragonite_saturation_state.values, **kwargs)  # summer omega

    # print('omega min/max')
    # print(f'gapfill {np.nanmin(ds_gapfill.aragonite_saturation_state.values)} / {np.nanmax(ds_gapfill.aragonite_saturation_state.values)}')
    # print(f'fall {np.nanmin(ds_fall_ph.aragonite_saturation_state.values)} / {np.nanmax(ds_fall_ph.aragonite_saturation_state.values)}')

    # summary_dict['gapfill_min'].append(np.round(np.nanmin(ds_gapfill.aragonite_saturation_state.values), 2))
    # summary_dict['gapfill_max'].append(np.round(np.nanmax(ds_gapfill.aragonite_saturation_state.values), 2))
    # summary_dict['fall_min'].append(np.round(np.nanmin(ds_fall_ph.aragonite_saturation_state.values), 2))
    # summary_dict['fall_max'].append(np.round(np.nanmax(ds_fall_ph.aragonite_saturation_state.values), 2))

    # format x-axis
    ax20.tick_params(axis='x', rotation=45)
    ax21.tick_params(axis='x', rotation=45)
    ax22.tick_params(axis='x', rotation=45)

    ax1.invert_yaxis()
    plt.subplots_adjust(left=0.06, right=0.93, bottom=0.1, top=0.9)

    sname = os.path.join(savedir, 'multipanel_xsection_2024-spring-summer.png')
    plt.savefig(sname, dpi=200)
    plt.close()

    # # write summary csv file
    # df = pd.DataFrame(summary_dict)
    # df.set_index('variable', inplace=True)
    # df.to_csv(os.path.join(savedir, 'temperature_chl_DO_pH_omega-2024-fall.csv'))


if __name__ == '__main__':
    savedir = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/multi_panel_xsection'
    main(savedir)
