#! /usr/bin/env python3

import numpy as np


def calc_mld_upper(df,ref_depth,deltaT,_z,_temp,_sal,lat,lon):
    #find the depth index closest to reference depth
    zdiff = _z - ref_depth
    dref=np.where((zdiff>0) & (zdiff<=1) & (~np.isnan(zdiff)))[0]
    
    # dref=np.argmin(np.abs(np.abs(_z)-ref_depth))
    if (dref.size!=0) and (~np.isnan(_z[dref[0]])):
        dref=dref[0]
    # if (np.abs(_z[dref]-ref_depth)<=1) and (~np.isnan(_z[dref])):
        
        pres=p_from_z(-_z,lat)
        SA = SA_from_SP(_sal,pres,lon,lat)
        #calculate intial rho
        ini_rho=pot_rho_t_exact(SA,_temp,pres,0)
        if (~np.isnan(SA[dref])) and (~np.isnan(_temp[dref])) and (~np.isnan(pres[dref])):
            # ini_rho=pot_rho_t_exact(SA[dref:],_temp[dref:],pres[dref:],0)
            #calculate rho profile using reference depth
            refdepth_rho=pot_rho_t_exact(SA[dref],_temp[dref],pres[dref],0)
            #calculate rho -deltaT reference depth
            refdepth_rho_deltaT =pot_rho_t_exact(SA[dref],_temp[dref]-deltaT,pres[dref],0)
            #calculate delta sigma
            delta_sig=np.abs(refdepth_rho_deltaT - refdepth_rho)
            #calculate difference between rho profile and rho at reference depth
            dens_diff=np.abs(ini_rho[dref:] - refdepth_rho)
            #grab first index where this is true
            if np.where(dens_diff>=delta_sig)[0].size>0:
                
                mld_idx=np.where(dens_diff>=delta_sig)[0][0]
                mldU=_z[dref:][mld_idx]
                qcU=1
                
            else:
                mldU=df.depth.iloc[-1]
                qcU=2 #profile doesn't surpass delta_sigma threshold
        else:
            mldU=np.nan
            qcU= 5 #one of T S P being nan cuased related density vars to be nan       
    else:
        mldU=np.nan
        qcU= 3 #profile doesn't have depth close enough to ref_depth or depth[ref_depth] is nan
    
    return mldU,qcU


def gap(prange):
    """
    :param prange: pressure range of the profile
    :return: the maximum allowable data gap for the profile
    """
    if prange < 20:
        gap_threshold = 8
    elif 20 <= prange < 50:
        gap_threshold = 10
    elif 50 <= prange < 200:
        gap_threshold = 25
    elif 200 <= prange < 500:
        gap_threshold = 50
    else:
        gap_threshold = 75

    return gap_threshold


def profile_mld_2way(dataframe, ref_depth=4, deltaT=0.5):
    """
    based off of de Boyer Montegut (2007) and Rudzin (2017)
    dataframe: single glider profile
    ref_depth: reference depth to calculate potential desnity at, default is 4m
    deltaT: change in temperature threshold, default is 0.5
    qc: 1 - profile is good and has mld
        2 - profile doesn't surpass delta_signma threshold
        3 - profile doesn't have a depth close enough to ref_depth
        4 - profile is less than 10m long
        5 - T, S, P being nan cuased related density vars to be nan
        6 - mldL < mldU, relatively unstratified
    """
    _z=dataframe.depth_interpolated.values
    _temp=dataframe.temperature.values
    _sal=dataframe.salinity.values
    lat=dataframe.lat.values
    lon=dataframe.lon.values

    if (~np.isnan(_sal).all()) and (~np.isnan(_temp).all()) and (~np.isnan(_z).all()):
        mldU,qcU = calc_mld_upper(dataframe,ref_depth,deltaT,_z,_temp,_sal,lat,lon)
        mldL,qcL = calc_mld_lower(dataframe,deltaT,_z,_temp,_sal,lat,lon)
        if mldL < mldU:
            qcU=6
            qcL=6
            mldU=dataframe.depth.iloc[-1]
            mldL=dataframe.depth.iloc[-1]
    else:
        mldU=np.nan
        mldL=np.nan
        qcU=5 # one of T S P being nan caused related density vars to be nan
        qcL=5

    return mldU, mldL, qcU, qcL


def profile_mld_n2(df, mld_var='density', zvar='pressure', qi_threshold=0.5):
    """
    Written by Sam Coakley and Lori Garzio, Jan 2022
    Calculates the Mixed Layer Depth (MLD) for a single profile as the depth of max Brunt‐Vaisala frequency squared
    (N**2) from Carvalho et al 2016 (https://doi.org/10.1002/2016GL071205). Calculates a Quality Index to determine
    if the water column is well-mixed, hence calculating MLD is not appropriate and the function will return
    MLD = np.nan, from Lorbacher et al, 2006 doi:10.1029/2003JC002157. When "QI > 0.8 a well defined MLD results.
    For QI in the range 0.5– 0.8, increased uncertainty of the profile interpretation becomes evident and with
    QI < 0.5 no mixed layer interpretation is possible." (from Section 3.4)
    :param df: depth profile in the form of a pandas dataframe
    :param mld_var: the name of the variable for which MLD is calculated, default is 'density'
    :param zvar: the name of the depth variable in the dataframe, default is 'pressure'
    :param qi_threshold: quality index threshold for determining well-mixed water, default is 0.5
    :return: the depth of the mixed layer in the units of zvar and the max buoyancy frequency in units of s-2
    """
    df.dropna(subset=[mld_var], inplace=True)
    pN2 = np.sqrt(9.81 / np.nanmean(df[mld_var]) * np.diff(df[mld_var], prepend=np.nan) / np.diff(df[zvar], prepend=np.nan)) ** 2
    if len(df) < 5:  # if there are <5 data bins, don't calculate MLD
        mld = np.nan
        maxN2 = np.nan
        qi = np.nan
    elif np.sum(~np.isnan(pN2)) < 3:  # if there are <3 values calculated for pN2, don't calculate MLD
        mld = np.nan
        maxN2 = np.nan
        qi = np.nan
    elif np.nanmax(np.diff(df[zvar])) > gap(np.nanmax(df[zvar]) - np.nanmin(df[zvar])):
        # if there is a gap in the profile that exceeds the defined threshold, don't calculate MLD
        mld = np.nan
        maxN2 = np.nan
        qi = np.nan
    else:
        pressure_range = [np.nanmin(df[zvar]), np.nanmax(df[zvar])]

        maxN2 = np.nanmax(pN2)
        mld_idx = np.where(pN2 == maxN2)[0][0]

        # if the code finds the first or last data point to be the max pN2, return nan
        if np.logical_or(mld_idx == 0, mld_idx == len(df) - 1):
            mld = np.nan
            maxN2 = np.nan
            qi = np.nan
        else:
            mld = np.nanmean([df[zvar][mld_idx], df[zvar][mld_idx + 1]])

            if mld < 5:
                # if MLD is <5, return nan
                mld = np.nan
                maxN2 = np.nan
                qi = np.nan
            elif np.logical_or(mld < pressure_range[0] + 2, mld > pressure_range[1] - 2):
                # if MLD is within 2 dbar of the top or bottom of the profile, return nan
                mld = np.nan
                maxN2 = np.nan
                qi = np.nan
            else:
                if qi_threshold:
                    # find MLD  1.5
                    mld15 = mld * 1.5
                    mld15_idx = np.argmin(np.abs(df[zvar] - mld15))

                    # Calculate Quality index (QI) from Lorbacher et al, 2006 doi:10.1029/2003JC002157
                    surface_mld_values = df[mld_var][0:mld_idx]  # values from the surface to MLD
                    surface_mld15_values = df[mld_var][0:mld15_idx]  # values from the surface to MLD * 1.5

                    qi = 1 - (np.std(surface_mld_values - np.nanmean(surface_mld_values)) /
                              np.std(surface_mld15_values - np.nanmean(surface_mld15_values)))

                    if qi < qi_threshold:
                        # if the Quality Index is < the threshold, this indicates well-mixed water so don't return MLD
                        mld = np.nan
                        maxN2 = np.nan

    return mld, maxN2, qi
