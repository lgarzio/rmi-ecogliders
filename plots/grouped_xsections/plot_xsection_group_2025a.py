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

    spring_ph = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2025/ru39-20250423T1535/ncei_pH/ru39-20250423T1535-delayed.nc')
    spring_do = xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2025/ru43-20250423T1533/ncei_dmon/ru43-20250423T1533-delayed.nc')
    earlysummer =  xr.open_dataset('/Users/garzio/Documents/rucool/Saba/gliderdata/2025/ru43-20250605T1606/ncei_dmon/ru43-20250605T1606-delayed.nc') 

    # make summary dataframe to export as csv
    summary_dict = dict(variable=['temperature', 'salinity', 'chl', 'DO_mgL', 'pH', 'omega'],
                        spring_min=[],
                        spring_max=[],
                        gapfill_min=[],
                        gapfill_max=[])

    # plot cross sections of all deployments: temperature, chl, DO, pH, omega
    fig, axs = plt.subplots(6, 2, figsize=(16, 18), sharex='col', sharey=True)
    plt.subplots_adjust(top=.96, bottom=0.08, right=.95, left=.06, hspace=0.05, wspace=0.05)
    ax1 = axs[0, 0]  # spring temperature
    ax2 = axs[0, 1]  # early summer temperature
    ax3 = axs[1, 0]  # spring salinity
    ax4 = axs[1, 1]  # early summer salinity
    ax5 = axs[2, 0]  # spring chl
    ax6 = axs[2, 1]  # early summer chl
    ax7 = axs[3, 0]  # spring DO
    ax8 = axs[3, 1]  # early summer DO
    ax9 = axs[4, 0]  # spring pH
    ax10 = axs[4, 1]  # early summer pH
    ax11 = axs[5, 0]  # spring omega
    ax12 = axs[5, 1]   # early summer omega

    # plot temperature
    kwargs = dict()
    kwargs['clabel'] = 'Temperature ({})'.format(r'$\rm ^oC$')
    kwargs['cmap'] = cmo.cm.thermal
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [7, 23]
    kwargs['title'] = 'Spring 2025'
    kwargs['ylims'] = [-4, 105]

    pf.xsection(fig, ax1, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.temperature.values, **kwargs)

    kwargs['title'] = 'Early Summer 2025'
    kwargs['ylabel'] = None
    kwargs['colorbar'] = True
    pf.xsection(fig, ax2, earlysummer.time.values, earlysummer.depth_interpolated.values, earlysummer.temperature.values, **kwargs)

    print('temperature min/max')
    print(f'spring {np.nanmin(spring_ph.temperature.values)} / {np.nanmax(spring_ph.temperature.values)}')
    print(f'gapfill {np.nanmin(earlysummer.temperature.values)} / {np.nanmax(earlysummer.temperature.values)}')
    summary_dict['spring_min'].append(np.round(np.nanmin(spring_ph.temperature.values), 2))
    summary_dict['spring_max'].append(np.round(np.nanmax(spring_ph.temperature.values), 2))
    summary_dict['gapfill_min'].append(np.round(np.nanmin(earlysummer.temperature.values), 2))
    summary_dict['gapfill_max'].append(np.round(np.nanmax(earlysummer.temperature.values), 2))

    # plot salinity
    kwargs = dict()
    kwargs['clabel'] = 'Salinity (PSU)'
    kwargs['cmap'] = cmo.cm.haline
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [29, 35]

    pf.xsection(fig, ax3, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.salinity.values, **kwargs)  # latesummer chl

    kwargs['ylabel'] = None
    kwargs['colorbar'] = True
    pf.xsection(fig, ax4, earlysummer.time.values, earlysummer.depth_interpolated.values, earlysummer.salinity.values, **kwargs)  # fall chl
    
    print('salinity min/max')
    print(f'spring {np.nanmin(spring_ph.salinity.values)} / {np.nanmax(spring_ph.salinity.values)}')
    print(f'gapfill {np.nanmin(earlysummer.salinity.values)} / {np.nanmax(earlysummer.salinity.values)}')
    summary_dict['spring_min'].append(np.round(np.nanmin(spring_ph.salinity.values), 2))
    summary_dict['spring_max'].append(np.round(np.nanmax(spring_ph.salinity.values), 2))
    summary_dict['gapfill_min'].append(np.round(np.nanmin(earlysummer.salinity.values), 2))
    summary_dict['gapfill_max'].append(np.round(np.nanmax(earlysummer.salinity.values), 2))

    # plot chl
    kwargs = dict()
    kwargs['clabel'] = 'Chlorophyll a (ug/L)'
    kwargs['cmap'] = cmo.cm.algae
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [2, 10]

    pf.xsection(fig, ax5, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.chlorophyll_a.values, **kwargs)  # latesummer chl

    kwargs['ylabel'] = None
    kwargs['colorbar'] = True
    pf.xsection(fig, ax6, earlysummer.time.values, earlysummer.depth_interpolated.values, earlysummer.chlorophyll_a.values, **kwargs)  # fall chl

    print('chla min/max')
    print(f'spring {np.nanmin(spring_ph.chlorophyll_a.values)} / {np.nanmax(spring_ph.chlorophyll_a.values)}')
    print(f'gapfill {np.nanmin(earlysummer.chlorophyll_a.values)} / {np.nanmax(earlysummer.chlorophyll_a.values)}')
    summary_dict['spring_min'].append(np.round(np.nanmin(spring_ph.chlorophyll_a.values), 2))
    summary_dict['spring_max'].append(np.round(np.nanmax(spring_ph.chlorophyll_a.values), 2))
    summary_dict['gapfill_min'].append(np.round(np.nanmin(earlysummer.chlorophyll_a.values), 2))
    summary_dict['gapfill_max'].append(np.round(np.nanmax(earlysummer.chlorophyll_a.values), 2))

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
    kwargs['xticklabels'] = None
    kwargs['xlabel'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [2, 9]
    kwargs['date_fmt'] = '%m-%d'

    pf.xsection(fig, ax7, spring_do.time.values, spring_do.depth_interpolated.values, spring_do.oxygen_concentration_shifted_mgL.values, **kwargs)  # latesummer DO

    kwargs['ylabel'] = None
    kwargs['colorbar'] = True
    pf.xsection(fig, ax8, earlysummer.time.values, earlysummer.depth_interpolated.values, earlysummer.oxygen_concentration_shifted_mgL.values, **kwargs)  # fall DO


    print('DO min/max')
    print(f'spring {np.nanmin(spring_do.oxygen_concentration_shifted_mgL.values)} / {np.nanmax(spring_do.oxygen_concentration_shifted_mgL.values)}')
    print(f'gapfill {np.nanmin(earlysummer.oxygen_concentration_shifted_mgL.values)} / {np.nanmax(earlysummer.oxygen_concentration_shifted_mgL.values)}')
    summary_dict['spring_min'].append(np.round(np.nanmin(spring_do.oxygen_concentration_shifted_mgL.values), 2))
    summary_dict['spring_max'].append(np.round(np.nanmax(spring_do.oxygen_concentration_shifted_mgL.values), 2))
    summary_dict['gapfill_min'].append(np.round(np.nanmin(earlysummer.oxygen_concentration_shifted_mgL.values), 2))
    summary_dict['gapfill_max'].append(np.round(np.nanmax(earlysummer.oxygen_concentration_shifted_mgL.values), 2))

    # pH
    kwargs = dict()
    kwargs['clabel'] = 'pH'
    kwargs['title'] = None
    kwargs['cmap'] = cmo.cm.matter
    kwargs['xlabel'] = None
    kwargs['xticklabels'] = None
    kwargs['colorbar'] = None
    kwargs['vlims'] = [7.6, 8.1]

    pf.xsection(fig, ax9, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.pH.values, **kwargs)  # latesummer pH

    kwargs['ylabel'] = None
    kwargs['colorbar'] = True
    placeholder = np.full(np.shape(earlysummer.time.values), np.nan)
    pf.xsection(fig, ax10, earlysummer.time.values, earlysummer.depth_interpolated.values, placeholder, **kwargs)  # fall pH

    print('pH min/max')
    print(f'spring: {np.nanmin(spring_ph.pH.values)} / {np.nanmax(spring_ph.pH.values)}')
    print(f'gapfill: none')

    summary_dict['spring_min'].append(np.round(np.nanmin(spring_ph.pH.values), 2))
    summary_dict['spring_max'].append(np.round(np.nanmax(spring_ph.pH.values), 2))
    summary_dict['gapfill_min'].append('nan')
    summary_dict['gapfill_max'].append('nan')

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

    pf.xsection(fig, ax11, spring_ph.time.values, spring_ph.depth_interpolated.values, spring_ph.aragonite_saturation_state.values, **kwargs)  # latesummer omega

    kwargs['ylabel'] = None
    kwargs['colorbar'] = True
    pf.xsection(fig, ax12, earlysummer.time.values, earlysummer.depth_interpolated.values, placeholder, **kwargs) 

    print('omega min/max')
    print(f'spring: {np.nanmin(spring_ph.aragonite_saturation_state.values)} / {np.nanmax(spring_ph.aragonite_saturation_state.values)}')
    print(f'gapfill: none')

    summary_dict['spring_min'].append(np.round(np.nanmin(spring_ph.aragonite_saturation_state.values), 2))
    summary_dict['spring_max'].append(np.round(np.nanmax(spring_ph.aragonite_saturation_state.values), 2))
    summary_dict['gapfill_min'].append('nan')
    summary_dict['gapfill_max'].append('nan')

    # format x-axis
    ax11.tick_params(axis='x', rotation=45)
    ax12.tick_params(axis='x', rotation=45)

    #ax1.invert_yaxis()
    plt.subplots_adjust(left=0.06, right=0.93, bottom=0.1, top=0.9)

    sname = os.path.join(savedir, 'multipanel_xsection_2025-spring-earlysummer.png')
    plt.savefig(sname, dpi=200)
    plt.close()

    # write summary csv file
    df = pd.DataFrame(summary_dict)
    df.set_index('variable', inplace=True)
    df.to_csv(os.path.join(savedir, 'gliderdata_summary-2025-spring.csv'))


if __name__ == '__main__':
    savedir = '/Users/garzio/Documents/rucool/Saba/RMI/reports/202507-progressreport'
    main(savedir)
