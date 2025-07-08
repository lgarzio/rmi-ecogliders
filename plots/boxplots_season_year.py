#!/usr/bin/env python

"""
Author: Lori Garzio on 7/2/2025
Last modified: 7/7/2025
Boxplots of surface- and bottom-averaged glider data for inshore, midshelf, and offshore.
The box limits extend from the lower to upper quartiles (25%, 75%), with a line at the median and a diamond symbol at
the mean. The whiskers extend from the box by 1.5x the inter-quartile range (IQR). Circles indicate outliers.
Notch indicates 95% CI around the median.
Statistics: Kruskal‐Wallis analysis with Dunn post hoc
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import yaml
import matplotlib.pyplot as plt
from scipy.stats import kruskal
import scikit_posthocs as sp
plt.rcParams.update({'font.size': 12})
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def build_summary(data, summary_dict, summary_dict_key):
    """
    :param data: data in the form of a numpy array
    :param summary_dict: dictionary to which summary data are appended
    :param summary_dict_key: key for dictionary
    :return:
    """
    if len(data) == 0:
        summary_dict[summary_dict_key] = dict(count=len(data),
                                              median=np.nan,
                                              mean=np.nan,
                                              lower_quartile=np.nan,
                                              upper_quartile=np.nan,
                                              lower_whisker=np.nan,
                                              upper_whisker=np.nan,
                                              min=np.nan,
                                              max=np.nan)
    else:
        lq = np.percentile(data, 25)
        uq = np.percentile(data, 75)
        iqr = uq - lq
        summary_dict[summary_dict_key] = dict(count=int(len(data)),
                                              median=np.round(np.nanmedian(data), 4),
                                              mean=np.round(np.nanmean(data), 4),
                                              lower_quartile=np.round(lq, 4),
                                              upper_quartile=np.round(uq, 4),
                                              lower_whisker=data[data >= lq - 1.5 * iqr].min(),
                                              upper_whisker=data[data <= uq + 1.5 * iqr].max(),
                                              min=np.round(np.nanmin(data), 4),
                                              max=np.round(np.nanmax(data), 4))


def get_variable_data(dataset, variable, loc):
    try:
        dataset1 = dataset[variable].sel(shelf_location=loc).values
        dataset1 = dataset1[~np.isnan(dataset1)]  # remove nans
    except KeyError:
        dataset1 = np.array([])
    return dataset1


def set_box_colors(bp, colors):
    for key in ['boxes', 'medians', 'fliers', 'means']:
        for patch, color in zip(bp[key], colors):
            patch.set_color(color)
            if key == 'boxes':
                patch.set_facecolor('none')
            elif key == 'means':
                patch.set_markerfacecolor(color)
                patch.set_markeredgecolor(color)
            elif key == 'fliers':
                patch.set_markeredgecolor(color)

    wc_colors = [x for pair in zip(colors, colors) for x in pair]
    for key in ['whiskers', 'caps']:
        for patch, color in zip(bp[key], wc_colors):
            patch.set_color(color)


def main(fname, save_dir):
    ds = xr.open_dataset(fname)
    ds = ds.swap_dims({'time': 'shelf_location'})
    surf_bot = fname.split('/')[-1].split('_')[2]

    save_dir = os.path.join(save_dir, surf_bot)
    os.makedirs(save_dir, exist_ok=True)
    summary_savedir = os.path.join(save_dir, 'summary_csv')
    os.makedirs(os.path.join(save_dir, 'summary_csv'), exist_ok=True)

    years = ['2023', '2024', '2025']
    box_colors = ['tab:blue', 'tab:orange', 'k']
    #box_colors = ['#440154', '#2A788E', '#7AD151']  # purple, blue, green

    # grab science vars to plot from config file
    root_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(root_dir)
    configdir = os.path.join(parent_dir, 'config')

    with open(os.path.join(configdir, 'boxplot_vars.yml')) as f:
        sci_vars = yaml.safe_load(f)

    # initialize empty dictionary for summary to export
    summary_dict = dict()

    # for each variable make boxplots and summary statistics
    for season in ['Winter', 'Spring', 'Summer', 'Fall']:
        season_idx = np.where(ds.season_year.str.contains(season))[0]
        ds_season = ds.isel(shelf_location=season_idx)
        unique_seasonyears = np.unique(ds_season.season_year)
        yr1_idx = np.where(ds_season.season_year.str.contains('2023'))[0]
        yr2_idx = np.where(ds_season.season_year.str.contains('2024'))[0]
        yr3_idx = np.where(ds_season.season_year.str.contains('2025'))[0]
        ds_season_yr1 = ds_season.isel(shelf_location=yr1_idx)
        ds_season_yr2 = ds_season.isel(shelf_location=yr2_idx)
        ds_season_yr3 = ds_season.isel(shelf_location=yr3_idx)
        for sv, info in sci_vars.items():
            inyr1 = get_variable_data(ds_season_yr1, sv, 'inshore')
            midyr1 = get_variable_data(ds_season_yr1, sv, 'midshelf')
            offyr1 = get_variable_data(ds_season_yr1, sv, 'offshore')
            inyr2 = get_variable_data(ds_season_yr2, sv, 'inshore')
            midyr2 = get_variable_data(ds_season_yr2, sv, 'midshelf')
            offyr2 = get_variable_data(ds_season_yr2, sv, 'offshore')
            inyr3 = get_variable_data(ds_season_yr3, sv, 'inshore')
            midyr3 = get_variable_data(ds_season_yr3, sv, 'midshelf')
            offyr3 = get_variable_data(ds_season_yr3, sv, 'offshore')

            # perform Kruskal-Wallis test for each season with Dunn post hoc
            labels = ['inshore_2023', 'midshelf_2023', 'offshore_2023',
                      'inshore_2024', 'midshelf_2024', 'offshore_2024',
                      'inshore_2025', 'midshelf_2025', 'offshore_2025']
            kw_list = [inyr1, midyr1, offyr1, inyr2, midyr2, offyr2, inyr3, midyr3, offyr3]
            kw_data = []
            kw_labels = []
            sample_sizes = []
            for i, d in enumerate(kw_list):
                if len(d) > 0:
                    kw_labels.append(labels[i])
                    kw_data.append(d)
                    sample_sizes.append(len(d))

            kw_stat, kw_p = kruskal(*kw_data, nan_policy='omit')
            if kw_p < 0.05:
                dunn = sp.posthoc_dunn(kw_data, p_adjust='bonferroni')
                dunn.columns = kw_labels
                dunn.index = kw_labels

                # add sample sizes
                dunn.loc['sample_size'] = sample_sizes

                stats_filename = f'{season}_{surf_bot}_{sv}_dunn_posthoc.csv'
                stats_savefile = os.path.join(summary_savedir, stats_filename)
                dunn.to_csv(stats_savefile)

            # append summary data to dictionary
            summary_dict[sv] = dict()
            build_summary(inyr1, summary_dict[sv], f'inshore_2023')
            build_summary(midyr1, summary_dict[sv], f'midshelf_2023')
            build_summary(offyr1, summary_dict[sv], f'offshore_2023')
            build_summary(inyr2, summary_dict[sv], f'inshore_2024')
            build_summary(midyr2, summary_dict[sv], f'midshelf_2024')
            build_summary(offyr2, summary_dict[sv], f'offshore_2024')
            build_summary(inyr3, summary_dict[sv], f'inshore_2025')
            build_summary(midyr3, summary_dict[sv], f'midshelf_2025')
            build_summary(offyr3, summary_dict[sv], f'offshore_2025')

            # make boxplot
            inshore_data = [list(inyr1), list(inyr2), list(inyr3)]
            midshelf_data = [list(midyr1), list(midyr2), list(midyr3)]
            offshore_data = [list(offyr1), list(offyr2), list(offyr3)]

            fig, ax = plt.subplots(figsize=(8, 7))

            # customize the boxplot elements
            meanpointprops = dict(marker='D')

            bp_in = ax.boxplot(inshore_data, positions=[1, 2, 3], widths=0.6, patch_artist=True, showmeans=True,
                                notch=True, meanprops=meanpointprops, sym='.')
            bp_mid = ax.boxplot(midshelf_data, positions=[5, 6, 7], widths=0.6, patch_artist=True, showmeans=True,
                                notch=True, meanprops=meanpointprops, sym='.')
            bp_off = ax.boxplot(offshore_data, positions=[9, 10, 11], widths=0.6, patch_artist=True, showmeans=True,
                                notch=True, meanprops=meanpointprops, sym='.')

            # set box colors
            set_box_colors(bp_in, box_colors)
            set_box_colors(bp_mid, box_colors)
            set_box_colors(bp_off, box_colors)

            # draw temporary lines and use them to create a legend
            plt.plot([], c='tab:blue', label='2023')
            plt.plot([], c='tab:orange', label='2024')
            plt.plot([], c='k', label='2025')
            plt.legend(fontsize=12)

            # set axes labels
            # ax.set_xticks([2, 6])
            # ax.set_xticklabels(['Surface', 'Bottom'])
            ax.set_xticks([2, 6, 10])
            ax.set_xticklabels(['Inshore', 'Midshelf', 'Offshore'])
            ax.set_ylabel(f'{ds[sv].attrs['long_name']} ({ds[sv].attrs['units']})')

            try:
                ylims = [info[surf_bot]['ymin'], info[surf_bot]['ymax']]
            except KeyError:
                ylims = ax.get_ylim()
            ax.set_ylim(ylims)
            ax.vlines(4, ylims[0], ylims[1], colors='dimgray', alpha=.5)
            ax.vlines(8, ylims[0], ylims[1], colors='dimgray', alpha=.5)

            ax.set_xlim([0, 12])
            ax.set_title(f'{season}: {surf_bot} {ds[sv].attrs['long_name']}')
            sfilename = f'boxplot_{sv}_{season}_{surf_bot}.png'
            sfile = os.path.join(save_dir, sfilename)
            plt.savefig(sfile, dpi=300)
            plt.close()

        # export summary as .csv
        df_list = []
        for k, v in summary_dict.items():
            dfk = pd.DataFrame(v)
            dfk.reset_index(inplace=True)
            dfk.index = list(np.repeat(k, len(dfk)))
            df_list.append(dfk)

        df = pd.concat(df_list, axis=0)
        df = df.round(2)
        df = df.rename(columns={'index': 'statistic'})
        csv_filename = f'{season}_{surf_bot}_boxplot_summary.csv'
        csv_savefile = os.path.join(summary_savedir, csv_filename)
        df.to_csv(csv_savefile)


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/data/RMI_gliderdata_surface_2023_2025.nc'
    save_directory = '/Users/garzio/Documents/rucool/Saba/RMI/plots_for_2025_report/boxplots'
    main(ncfile, save_directory)
