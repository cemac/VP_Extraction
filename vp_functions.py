# CVP functional code
# Used to extract CVP data from cf-radial files using Pyart
# Adpated by R. Neely March 2022

import numpy as np
import pyart
from netCDF4 import num2date
import kdp_functions as kdpfun
import math
import datetime as dt
from scipy import stats

from vp import *


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

'''
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except (ValueError):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
'''

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.


    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    box = np.ones(window_len)/window_len
    y = np.convolve(x, box, mode='same')
    y[0] = x[0]# to remove jumps at the lowest elevation
    return y


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
def altitude_parameter_averaging_cvp(radar, vp, field, tix, timeofsweep, rixs, azixs, all_bin_indexes, azimuth_exclude = None):

    field_fill_to_nan(radar, field)
    nazimuths=int(radar.nrays/radar.nsweeps)
    data = np.array(radar.fields[field]['data'].reshape(radar.nsweeps,nazimuths,radar.ngates))

    column = data[np.ix_(range(0,radar.nsweeps),azixs,rixs)]


    if field in ['dBuZ', 'dBZ', 'dBZ_ac', 'dBuZv', 'dBZv']:column = np.power(10, column / 10.0)

    # unfold uPhiDP values
    # remove phi_dp wrap-around:
    #print 'unwrapping phidp ...'
    if field in ['uPhiDP']:
        elev = column.shape[0]
        rays = column.shape[1]
        bins = column.shape[2]
        METEO_THRESH=0.7
        flags = np.zeros((elev, rays, bins))
        # generate non-meteo mask:
        rhohv_data = radar.fields['RhoHV']['data']
        rhohv_3D_data = rhohv_data.reshape((radar.nsweeps,nazimuths,radar.ngates))
        rhohv = rhohv_3D_data[np.ix_(range(0,radar.nsweeps),azixs,rixs)]
        (meteoMask) = kdpfun.generate_meteo_mask(elev, rays, bins, flags, rhohv, METEO_THRESH)
        column = kdpfun.unwrap_phidp(elev, rays, bins, meteoMask, column)

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

    if len(vp.heights) >= 300: win = 5
    else:win = 3

    equdist_mean = smooth(equdist_mean,win)

    if field in ['dBuZ', 'dBZ', 'dBZ_ac', 'dBuZv', 'dBZv']:
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
                                     meteoMask=[], snrMask=[], verbose=False):

    field_fill_to_nan(radar, field)
    nazimuths=int(radar.nrays/radar.nsweeps)
    data = np.array(radar.fields[field]['data'].reshape(radar.nsweeps,nazimuths,radar.ngates))
    sweep_data=data[elev_ix, :,:] # leave as 3d but just one elevation
    if len(meteoMask)>0:
        sweep_data[(meteoMask != 1) | (snrMask != 1)] = np.nan

    # unfold uPhiDP values
    # remove phi_dp wrap-around:
    if field in ['uPhiDP']:
        sweep_data = kdpfun.unwrap_phidp(1, nazimuths, radar.ngates, meteoMask, sweep_data)
    elif field in ['ZDR'] and len(meteoMask)>0:
        sweep_data[(sweep_data < -1.5) | (sweep_data > 5)] = np.nan
    elif field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:
        if len(meteoMask)>0:
            sweep_data[(sweep_data < -10) | (sweep_data > 60)] = np.nan
        sweep_data = np.power(10, sweep_data / 10.0)

    sweep_data = np.ma.masked_where(azimuth_mask,sweep_data)

    mean_values=np.nanmean(sweep_data, axis=1)[0,:]
    std_values=np.nanstd(sweep_data, axis=1)[0,:]
    counts=np.sum(np.isfinite(sweep_data), axis=1).data[0,:]
    if field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:
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
                    meteo=False,
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
    
    if meteo:
        flags = np.zeros((radar.nsweeps,nazimuths,radar.ngates))
        rhohv = radar.fields['RhoHV']['data'].reshape((radar.nsweeps,nazimuths,radar.ngates))[elev_ix,:,:]
        METEO_THRESH=0.7
        (meteoMask) = kdpfun.generate_meteo_mask(1,nazimuths,radar.ngates, flags, rhohv, METEO_THRESH)

        snrH = radar.fields['SNR']['data'].reshape(radar.nsweeps,nazimuths,radar.ngates)[elev_ix,:,:]
        snrV = radar.fields['SNRv']['data'].reshape(radar.nsweeps,nazimuths,radar.ngates)[elev_ix,:,:]
        snrMask = np.full_like(meteoMask,0)
        snrMask[(snrH > 8) & (snrV > 8)] = 1
    else:
        meteoMask=[]
        snrMask=[]
    
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


