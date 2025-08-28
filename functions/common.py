#! /usr/bin/env python3

"""
Author: Lori Garzio on 6/24/2025
Last modified: 8/28/2025
"""

import numpy as np
import pandas as pd
from cartopy.io.shapereader import Reader
from shapely.ops import unary_union
from shapely.geometry import Polygon


def depth_bin(dataframe, depth_var='depth', depth_min=0, depth_max=None, stride=1):
    """
    Written by Mike Smith
    :param dataframe: depth profile in the form of a pandas dataframe
    :param depth_var: the name of the depth variable in the dataframe
    :param depth_min: the shallowest bin depth
    :param depth_max: the deepest bin depth
    :param stride: the amount of space between each bin
    :return: pandas dataframe where data has been averaged into specified depth bins
    """
    depth_max = depth_max or dataframe[depth_var].max()

    bins = np.arange(depth_min, depth_max+stride, stride)  # Generate array of depths you want to bin at
    cut = pd.cut(dataframe[depth_var], bins)  # Cut/Bin the dataframe based on the bins variable we just generated
    binned_df = dataframe.groupby(cut, observed=False).mean()  # Groupby the cut and do the mean
    return binned_df


def return_noaa_polygons(extend=False):
    """
    Combine multiple strata levels defined in the NOAA bottom trawl survey strata downloaded from
    https://github.com/NOAA-EDAB/FisheryConditionLinks/tree/master/NES_BOTTOM_TRAWL_STRATA
    for the NYB region into inshore, midshelf, and offshore. STRATA mapping from Laura Nazzaro.
    :param extend: if True, return extended polygons that include regions farther south
    """
    shpfile = Reader(
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/NES_BOTTOM_TRAWL_STRATA/NES_BOTTOM_TRAWL_STRATA.shp')
    shplist = list(shpfile.records())

    strata_mapping = dict(
        inshore=dict(
            snum=[3150, 3160, 3170, 3180, 3190, 3200, 3120, 3130, 3140, 3100, 3090, 3110, 3060, 3070, 3080],
            color='#9e36d7ff',
            zorder=20,
            poly=[]
        ),  # purple
        midshelf=dict(snum=[1730, 1010], color='#d7369eff', zorder=19, poly=[]),  # pink
        offshore=dict(snum=[1740, 1750, 1760, 1020, 1030, 1040], zorder=18, color='#95c983ff', poly=[])  # green
    )

    if extend:
        strata_mapping = dict(
            inshore=dict(
                #snum=[3150, 3160, 3170, 3180, 3190, 3200, 3230, 3120, 3130, 3140, 3100, 3090, 3110, 3060, 3070, 3080],
                snum=[3120, 3130, 3140, 3100, 3090, 3110, 3060, 3070, 3080, 3150, 3160, 3170, 3180, 3190, 3200, 3210, 
                      3220, 3230, 3240, 3250, 3260, 3270, 3280, 3290, 3300, 3310, 3320, 3030, 3050, 3040, 3010, 3020, 3450],
                color='#440154',  # purple
                zorder=20,
                poly=[]
            ), 
            midshelf=dict(
                #snum=[1730, 1010, 1690],
                snum=[1010, 1730, 1690, 1050],
                color='cyan',  # '#2A788E' blue
                zorder=19, 
                poly=[]
            ),  
            offshore=dict(
                #snum=[1740, 1750, 1760, 1020, 1030, 1040, 1700, 1710, 1720], 
                snum=[1020, 1030, 1040, 1740, 1750, 1760, 1700, 1710, 1720, 1060, 1070, 1080],
                color='#7AD151',  # green
                zorder=18, 
                poly=[]
            ) 
        )

    # combine regions to NYB inshore, midshelf, offshore
    for region_code, sm in strata_mapping.items():
        polys = []
        for sl in shplist:
            if sl.attributes['STRATA'] in sm['snum']:
                poly = sl.geometry
                polys.append(poly)
        if region_code == 'offshore':  # add a little extra polygon to the offshore region
            lon = [-72.85, -72.85, -73.2, -73.2, -72.85]
            lat = [38.86, 38.75, 38.75, 38.86, 38.86]
            extra_poly = Polygon(zip(lon, lat))
            polys.append(extra_poly)

        outside_poly = unary_union(polys)
        strata_mapping[region_code]['poly'] = outside_poly

    return strata_mapping
