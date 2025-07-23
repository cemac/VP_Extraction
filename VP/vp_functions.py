# CVP functional code
# Used to extract CVP data from cf-radial files using Pyart
# Adpated by R. Neely March 2022

import numpy as np
import pyart
from netCDF4 import num2date
import math
import datetime as dt
from scipy import stats

from .vp import VerticalProfile


def field_fill_to_nan(radar, mask_field):

    if mask_field in radar.fields.keys():
        try:
            if(np.ma.is_masked(radar.fields[mask_field]['data'])):
                mask = radar.fields[mask_field]['data'].mask
                radar.fields[mask_field]['data'] = radar.fields[mask_field]['data'].data
                radar.fields[mask_field]['data'][mask==True] = np.nan
        except:
            print('''field_fill_to_nan doesn't work, try again''')
            raise
    else:
        print ("{} is not in radar.fields.keys() {}".format(mask_field,radar.fields.keys()))
    return

    
#------------------------------------------------------------------------------------------------------
# get the centre rix and az ix in the radar data at elevation index elev_ix, for the point that most closely matches the lat and lon 
#------------------------------------------------------------------------------------------------------
def get_centre_azix_rix_for_lat_lon(radar, elev_ix, lat, lon):

    # don't assume that elevations are in order lowest to highest so find index of lowest
    nazimuths=int(radar.nrays/radar.nsweeps)

    lat_data = radar.gate_latitude['data'].reshape(radar.nsweeps,nazimuths,radar.ngates)
    lon_data = radar.gate_longitude['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))
    X = np.sqrt( np.square( lat_data[elev_ix,:,:] - lat ) +  np.square( lon_data[elev_ix,:,:] - lon ) )
    idx = np.where( X == X.min() )

    rix = idx[1][0]
    azix = idx[0][0]
    return rix, azix
    
#------------------------------------------------------------------------------------------------------
# get the lat and lon at the centre of the radar data 
#------------------------------------------------------------------------------------------------------
def get_centre_lat_lon_for_radar(radar):

    nazimuths=int(radar.nrays/radar.nsweeps)
    lat_data = radar.gate_latitude['data'].reshape(radar.nsweeps,nazimuths,radar.ngates)
    lon_data = radar.gate_longitude['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))
    centre_lat=np.mean(lat_data[:,:,0])
    centre_lon=np.mean(lon_data[:,:,0])

    return centre_lat, centre_lon


#----------------------------------------------------------------------------
# get all the azimuth indices for a CVP that is centred on centre_rix, centre_azix with radius col_radius
#----------------------------------------------------------------------------
def get_az_indexes(centre_rix, centre_azix, col_radius, radar, verbose = True):

    # delta azimuths are 360/naz and naz=radar.nrays/radar.nsweeps
    range_dist=radar.range['data'][centre_rix]-radar.range['data'][0]
    naz=int(radar.nrays/radar.nsweeps)
    delta_az=360.0/naz
    tan_half_beam=math.tan(math.radians(delta_az/2.0))
    az_size = 2.0 * range_dist * tan_half_beam
    avg_azix_delta = int(col_radius * 1000/az_size)
#    if(verbose):
#        print ("get_az_indexes(): azimuthal size for this location is {}".format(az_size))
#        print ("int({} * 1000/{}) = avg_azix_delta is {}".format(col_radius,az_size,avg_azix_delta))

    if(centre_azix<avg_azix_delta):
        from_azix = range(naz + (centre_azix-avg_azix_delta),naz)
        to_azix = range(0,centre_azix+avg_azix_delta+1)
    elif(centre_azix+avg_azix_delta>naz-1):
        from_azix = range(0,(centre_azix+avg_azix_delta)-(naz-1))
        to_azix = range(centre_azix-avg_azix_delta,naz)
    else:
        from_azix = range(centre_azix-avg_azix_delta,centre_azix)
        to_azix = range(centre_azix,centre_azix+avg_azix_delta+1)
    az_indxs = []
    az_indxs.extend(from_azix)
    az_indxs.extend(to_azix)
    return(az_indxs)


#----------------------------------------------------------------------------
# get all the range indices for a CVP that is centred on centre_rix, with radius col_radius
#----------------------------------------------------------------------------
def get_r_indexes(centre_rix, radar, col_radius, verbose = True):

    steps_in_range_delta = int(col_radius * 1000/(radar.range['data'][1]-radar.range['data'][0]))
    max_rix=len(radar.range['data'])
#    if (verbose):
#        print ("in get_r_indexes(centre_rix = {},steps_in_range_delta = {}, max rix= {})".format(centre_rix, steps_in_range_delta, max_rix))

    if(centre_rix+steps_in_range_delta>max_rix):
        from_rix = range(centre_rix-steps_in_range_delta-1,centre_rix)
        to_rix =  range(centre_rix,max_rix)
    elif(centre_rix-steps_in_range_delta<0):
        from_rix = range(0,centre_rix)
        to_rix = range(centre_rix,centre_rix+steps_in_range_delta)
    else:
        from_rix = range(centre_rix-steps_in_range_delta,centre_rix)
        to_rix = range(centre_rix,centre_rix+steps_in_range_delta+1)

    r_indxs = []
    r_indxs.extend(from_rix)
    r_indxs.extend(to_rix)
    return(r_indxs)

#----------------------------------------------------------------------------
# Perform the averaging over the column defined by vp (VerticalProfile object) for requested field
# This sets up the data in the vp object for time index tix
# For a static CVP we would expect the lat and lon and bin_indexes to be the same for each time whereas
# a dynamic cvp will have changing lat and lon and the altitudes in the different elevations may be
# different so bin_indexes can be different each time
# inputs:
#    radar - a radar object from which we get the data
#    field - name of field to get
#    tix -  the time index to use when setting up the data in the vp
#    timeofsweep - timestamp for this data
#    all_bin_indexes - list of bin indexes of length nheights - gives the elevation, azimuth and range indices
#                                                               for this column for this time
#    azimuth_exclude - this currently does not do anything with azimuth exclude
#----------------------------------------------------------------------------
def altitude_parameter_averaging_cvp(radar,
                                     vp,
                                     field,
                                     tix,
                                     timeofsweep,
                                     rixs,
                                     azixs,
                                     all_bin_indexes, 
                                     azimuth_exclude = None,
                                     logarithmic_field=False):

    field_fill_to_nan(radar, field)
    nazimuths=int(radar.nrays/radar.nsweeps)
    data = np.array(radar.fields[field]['data'].reshape(radar.nsweeps,nazimuths,radar.ngates))

    column = data[np.ix_(range(0,radar.nsweeps),azixs,rixs)]


    if logarithmic_field:
        column = np.power(10, column / 10.0)

    # get mean of means for equidistant
    equdist_mean = np.zeros(vp.heights.shape[0])*np.nan
    equdist_std = np.zeros(vp.heights.shape[0])*np.nan
    equdist_count = np.zeros(vp.heights.shape[0])
    equdist_total = np.zeros(vp.heights.shape[0])
    for h in range(vp.heights.shape[0]):
        bin_indexes=all_bin_indexes[h]
        if len(bin_indexes[0])>0:
            temp_column = column[bin_indexes]
            ix=np.where(np.isfinite(temp_column))
            if len(ix[0])>0:
                equdist_mean[h] = np.nanmean(temp_column[ix])
                equdist_std[h] = np.std(temp_column[ix])
                equdist_count[h] = len(ix[0])
                equdist_total[h] = len(temp_column)

    if logarithmic_field:
        equdist_mean = 10 * np.log10(equdist_mean)
        equdist_std = 10 * np.log10(equdist_std)


    radar_elevations = radar.elevation['data'].reshape((radar.nsweeps,nazimuths))
    radar_elevations=np.unique(radar_elevations)
    emin_ix=np.argmin(radar_elevations)    
    lat_data = radar.gate_latitude['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))
    lon_data = radar.gate_longitude['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))
    nr=len(rixs)
    naz=len(azixs)
    centre_rix=rixs[int(nr/2)]
    centre_azix=azixs[int(naz/2)]
    column_lats=lat_data[:,centre_azix,centre_rix]
    column_lons=lon_data[:,centre_azix,centre_rix]
    lat=column_lats[emin_ix]
    lon=column_lons[emin_ix]

    vp.add_data_for_time(field, tix, timeofsweep, equdist_mean, equdist_std, equdist_count, lon, lat)


def altitude_parameter_averaging_qvp(radar, vp, elev_ix, field, tix, timeofsweep, azimuth_mask,
                                     meteoMask=[], snrMask=[],
                                     logarithmic_field=False,
                                     verbose=False):

    field_fill_to_nan(radar, field)
    nazimuths=int(radar.nrays/radar.nsweeps)
    data = np.array(radar.fields[field]['data'].reshape(radar.nsweeps,nazimuths,radar.ngates))
    sweep_data=data[elev_ix, :,:] # leave as 3d but just one elevation
    if len(meteoMask)>0:
        sweep_data[(meteoMask != 1) | (snrMask != 1)] = np.nan

    
    if logarithmic_field:
        sweep_data = np.power(10, sweep_data / 10.0)

    sweep_data = np.ma.masked_where(azimuth_mask,sweep_data)

    mean_values=np.nanmean(sweep_data, axis=1)[0,:]
    std_values=np.nanstd(sweep_data, axis=1)[0,:]
    counts=np.sum(np.isfinite(sweep_data), axis=1).data[0,:]
    if logarithmic_field:
        mean_values = 10 * np.log10(mean_values)
        std_values = 10 * np.log10(std_values)

    # clean up the nearest beans influenced by sidelobs
    ranges=radar.range['data']
    mean_values[ranges<400] = np.nan

    lat_data = radar.gate_latitude['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))
    lon_data = radar.gate_longitude['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))
    # use closest range lat and lons to get centre values
    lat=np.mean(lat_data[elev_ix,:,0])
    lon=np.mean(lon_data[elev_ix,:,0])
    vp.add_data_for_time(field, tix, timeofsweep, mean_values, std_values, counts, lon, lat)

def time_height_qvp(radar, vp, field_list, tix,
                    elevation, azimuth_exclude = [],
                    meteoMask=None,
                    snrMask=None,
                    verbose=False):

    timeofsweep = num2date(np.nanmean(radar.time['data'][:]),
                                       radar.time['units'],
                                       radar.time['calendar'])
    nazimuths=int(radar.nrays/radar.nsweeps)
    elevations=radar.elevation['data'].reshape(radar.nsweeps, nazimuths)[:,0] # only need to look at one azimuth as they all have the same elevation
    elev_ix = np.where(abs(elevations - elevation) <= 0.1)[0]
    if len(elev_ix)!=1:
        raise ValueError('time_height_qvp: cannot find elevation {} in radar data'.format(elevation))
        
    # create azimuth exclude mask
    azimuth_data = radar.azimuth['data'].reshape(radar.nsweeps, nazimuths)[elev_ix,:]
    azimuth_mask = np.zeros((1,nazimuths, radar.ngates), bool)
    for bounds in azimuth_exclude:
        ix=np.where((azimuth_data>=bounds[0])&(azimuth_data<=bounds[1]))
        azimuth_mask[ix]=True
    inv_mask = np.where(azimuth_mask, 0, 1)
    max_counts=np.sum(inv_mask, axis=1)
    vp.max_counts[:, tix]=max_counts
    
    for field in field_list:
        altitude_parameter_averaging_qvp(radar, vp, elev_ix, field, tix, timeofsweep,
                                         azimuth_mask,
                                         meteoMask=meteoMask, snrMask=snrMask,
                                         verbose=verbose)
        
def time_height_cvp(radar, vp, col_radius, field_list, tix,
                    equidistant_bound,
                    azimuth_exclude = [],
                    verbose=False):

    timeofsweep = num2date(np.nanmean(radar.time['data'][:]),
                                       radar.time['units'],
                                       radar.time['calendar'])
    nazimuths=int(radar.nrays/radar.nsweeps)
    # get altitudes from the radar object
    altitudes = radar.fields['scan_altitude']['data']#[radar.sweep_start_ray_index['data'][sweep], :]

    # reshape the field to 3D array
    altitudes = altitudes.reshape((radar.nsweeps,nazimuths,radar.ngates))

    # find the column indices
    lat_data = radar.gate_latitude['data'].reshape(radar.nsweeps,nazimuths,radar.ngates)
    lon_data = radar.gate_longitude['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))

    # for now use the minimum elevation to get rix and azix - we should really do this for each elevation separately
    # or we should find any voxels  that are withon the lat lon bounds
    # but this is how it used to be
    radar_elevations = radar.elevation['data'].reshape((radar.nsweeps,nazimuths))
    radar_elevations=radar_elevations[:,0] # elevations are same for all azimuths
    eix=np.argmin(radar_elevations)
    rix, azix=get_centre_azix_rix_for_lat_lon(radar, eix, vp.centre_lat_lon[0], vp.centre_lat_lon[1])
    # get the sector indices around the centre indices
    azixs=get_az_indexes(rix, azix, col_radius, radar, verbose=verbose)
    rixs=get_r_indexes(rix, radar, col_radius, verbose = verbose)
    # get altitudes for the column values
    column_altitudes = altitudes[np.ix_(np.arange(radar.nsweeps),azixs,rixs)]

    # find the indices for each height
    all_bin_indexes=[]
    for item, boundary in enumerate(equidistant_bound[0:-1]):

        if(item == len(equidistant_bound)-1):
            bin_indexes = np.where(column_altitudes>=boundary)
        else:
            bin_indexes = np.where((column_altitudes>=boundary) & (column_altitudes<equidistant_bound[item+1]))
        all_bin_indexes.append(bin_indexes)
        vp.max_counts[item, tix]=len(bin_indexes[0])

    # now calculate the means, stddev etc. for each field
    for field in field_list:

        altitude_parameter_averaging_cvp(radar, vp, field, tix, timeofsweep,
                                         rixs, azixs, all_bin_indexes,
                                         azimuth_exclude = azimuth_exclude)
        
# ===========================================================================================================
# CODE FOR LOCATING COLUMNS - 2D masks returned - column_mask
# Current implementation - rectangular cartesian column
#                        - rectangular cartesian column from lat/lon
#                        - polar column (degrees in azi, metres in range) from an azi, horizontal range
#                        - polar column from lat/lon
#                        -
# TO ADD IN FUTURE - Cylindrical columns using a radius
# ===========================================================================================================

def generate_cartesian_column_mask(radar, target_x_position, target_y_position, x_size, y_size):
    """
    Return a 2D boolean mask that matches the dimensions of the radar data (nrays,ngates).
    The column will be centred on the target position with the x and y size being independently specified,
    which means rectangular columns could be extracted if desired.
    """
    column_mask = np.all([radar.gate_x['data']>target_x_position-x_size/2.0,
                                radar.gate_x['data']<target_x_position+x_size/2.0,
                                radar.gate_y['data']>target_y_position-y_size/2.0,
                                radar.gate_y['data']<target_y_position+y_size/2.0,
                                ],
                        axis=0)
    return column_mask

def generate_cartesian_column_mask_from_lat_lon(radar, target_longitude, target_latitude, x_size, y_size):
    """
    Return a 2D boolean mask that matches the dimensions of the radar data (nrays,ngates).
    The column will be centred on the target latitude and longitude with the x and y size being independently specified in cartesian units.
    """
    x,y = pyart.core.geographic_to_cartesian_aeqd(lon=target_longitude,lat=target_latitude,
                                        lon_0=radar.longitude['data'][0],
                                        lat_0=radar.latitude['data'][0],
                                       )
    
    column_mask = generate_cartesian_column_mask(radar, x, y, x_size, y_size)
    
    return column_mask
    
def _find_azimuth_mask(radar, target_azimuth, azi_step):
    """
    Return a 2D boolean mask containing radar azimuths within a given window (+/- the azimuth step) centred
    on the target azimuth.
    The mask will cover all sweeps within the radar object.    
    """
    azimuths = radar.azimuth['data']
    
    if target_azimuth+azi_step >=360.0:
        azi_mask = np.logical_or(azimuths>=target_azimuth-azi_step, 
                                 azimuths<=target_azimuth+azi_step-360)
    elif target_azimuth-azi_step < 0:
        azi_mask = np.logical_or(azimuths>=360+target_azimuth-azi_step,
                                 azimuths<=target_azimuth+azi_step)
    
    else:
        azi_mask = np.logical_and(azimuths>=target_azimuth-azi_step,
                                  azimuths<=target_azimuth+azi_step)

    azi_mask = np.repeat(azi_mask,radar.ngates).reshape(radar.nrays,radar.ngates)
    return(azi_mask)

def _find_range_mask(radar, target_range, range_step):

    horizontal_ranges = np.sqrt(radar.gate_x['data']**2 + radar.gate_y['data']**2)
    range_mask = np.logical_and(horizontal_ranges>=target_range-range_step,
                                horizontal_ranges<=target_range+range_step,
                                )
    
    return range_mask

def generate_polar_column_mask(radar, target_azimuth, target_range, azimuth_size, range_size):
    """
    Return a 2D boolean mask that matches the dimensions of the radar data (nrays,ngates).
    The column will be centred on the target azimuth and range with the azimuth and range size being independently specified (in degrees and metres),
    which means rectangular columns could be extracted if desired.
    """

    azimuth_mask = _find_azimuth_mask(radar, target_azimuth, azimuth_size/2.0)
    range_mask = _find_range_mask(radar, target_range, range_size/2.0)
    column_mask = np.all([azimuth_mask,
                          range_mask],
                        axis=0)
    return column_mask

def generate_polar_column_mask_from_lat_lon(radar, target_longitude, target_latitude, azimuth_size, range_size):
    """
    Return a 2D boolean mask that matches the dimensions of the radar data (nrays,ngates).
    The column will be centred on the target latitude and longitude with the azimuth and range size being independently specified (in degrees and metres),
    which means rectangular columns could be extracted if desired.
    """
    target_azimuth = pyart.util.for_azimuth(radar.latitude['data'],
                                              target_latitude,
                                              radar.longitude['data'],
                                              target_longitude)
    target_range = pyart.util.sphere_distance(radar.latitude['data'],
                                              target_latitude,
                                              radar.longitude['data'],
                                              target_longitude)
    azimuth_mask = _find_azimuth_mask(radar, target_azimuth, azimuth_size/2.0)
    range_mask = _find_range_mask(radar, target_range, range_size/2.0)
    column_mask = np.all([azimuth_mask,
                          range_mask],
                        axis=0)
    return column_mask

def combine_column_mask_and_gatefilter(column_mask, gatefilter):
    """
    Combine a CVP column mask with a pyart gatefilter object to apply filtering to the CVP column without 
    editing the underlying data within the radar object.
    """
    filtered_column_mask = np.logical_and(column_mask,
                                          gatefilter.gate_included)
    return filtered_column_mask