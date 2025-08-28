#! /usr/bin/env python

"""
Author: Lori Garzio on 2/14/2025
Last modified: 8/28/2025
"""

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams.update({'font.size': 12})


def map_add_shapefile_outlines(ax, shpfile, edgecolor=None, facecolor=None, projection=None, zorder=None, alpha=None,
                          linewidth=1):
    edgecolor = edgecolor or 'black'
    facecolor = facecolor or 'none'
    projection = projection or ccrs.PlateCarree()
    zorder = zorder or 1
    alpha = alpha or 1

    shape_feature = cfeature.ShapelyFeature(
        Reader(shpfile).geometries(),
        projection,
        edgecolor=edgecolor,
        facecolor=facecolor,
        linewidth=linewidth
        )

    ax.add_feature(shape_feature, zorder=zorder, alpha=alpha)

    return shape_feature._geoms


def xsection(fig, ax, x, y, z, xlabel='Time', xticklabels=True, ylabel='Depth (m)', clabel=None, cmap='jet', title=None, date_fmt=None,
             grid=None, extend='both', vlims=None, colorbar=True, ylims=None):
    if vlims:
        xc = ax.scatter(x, y, c=z, cmap=cmap, s=10, vmin=vlims[0], vmax=vlims[1], edgecolor='None')
    else:
        xc = ax.scatter(x, y, c=z, cmap=cmap, s=10, edgecolor='None')

    if not xticklabels:
        ax.tick_params(labelbottom=False)

    if ylims:
        ax.set_ylim(ylims)

    ax.invert_yaxis()

    if ylabel:
        ax.set_ylabel(ylabel)

    if xlabel:
        ax.set_xlabel(xlabel)

    if title:
        ax.set_title(title, fontsize=22)

    # format colorbar
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
        fig.add_axes(cax)
        if clabel:
            cb = plt.colorbar(xc, cax=cax, label=clabel, extend=extend)
        else:
            cb = plt.colorbar(xc, cax=cax, extend=extend)

    # format x-axis
    if date_fmt:
        xfmt = mdates.DateFormatter(date_fmt)
        ax.xaxis.set_major_formatter(xfmt)

    if grid:
        ax.grid(ls='--', lw=.5)
