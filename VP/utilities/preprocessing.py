# Functions for preprocessing QVP files, wrapped into a single preprocessing function
# Written by David Dufton, Nov. 2016

import numpy as np
import copy
import kdp_functions as kdpfun


#import nrt.corrections.attenuation as nrt_attenuation

FILTER_BINS_1 = 11      # median filter 1 window length
FILTER_BINS_2 = 9       # median filter 2 window length
METEO_THRESH = 0.7      # (rhohv) - remove non-meteo
RAIN_THRESH = 0.85      # (rhohv) - remove non-rain
SMOOTH_BINS_1 = 5       # smoothing filter 1 window length
SMOOTH_BINS_2 = 3       # smoothing filter 2 window length


def beam_height(r, e, h_):
    height = h_ + (r * np.sin(np.deg2rad(e))) + ((r ** 2) / (2 * (4 / 3.0) * 6371 * 1000))
    return height


def height_array(ranges, elevations, h_):
    h_arr = np.zeros((elevations.shape[0], ranges.shape[0]))
    for i, elevation in enumerate(elevations):
        h_arr[i] = beam_height(ranges, elevation, h_)
    return h_arr


def add_beam_height(radar):
    if radar.range['data'].ndim == 1:
        test = height_array(radar.range['data'], radar.elevation['data'], radar.altitude['data'])
    else:
        print ('Warning, radar range has more than one dimension')
        test = height_array(radar.range['data'][0], radar.elevation['data'], radar.altitude['data'])
    radar.add_field('scan_altitude', {'data': test, 'units': 'metres'}, replace_existing=True)

def kdp_ukmo(radar,
             phidpfield='uPhiDP',
             rhohvfield='RhoHV',
             FILTER_BINS_1=11,
             FILTER_BINS_2=9,
             METEO_THRESH=0.7,
             RAIN_THRESH=0.85,
             SMOOTH_BINS_1=5,
             SMOOTH_BINS_2=3):

    elev = 1
    rays = copy.deepcopy(radar.nrays)
    bins = copy.deepcopy(radar.ngates)
    phidp = copy.deepcopy(radar.fields[phidpfield]['data']).reshape(elev,rays,bins)
    rhohv = copy.deepcopy(radar.fields[rhohvfield]['data']).reshape(elev,rays,bins)
    binlength = radar.range['data'][1] - radar.range['data'][0]
    
    #fields stored as 2d arrays in radar object - don't need column code for cvp

    # flags = np.where(rad2.fields['classification']['data']==1,0,1)
    flags = np.zeros((elev, rays, bins))

    # generate non-meteo mask:
    #print 'generating non-meteo mask ...'
    (meteoMask) = kdpfun.generate_meteo_mask(elev, rays, bins, flags, rhohv, METEO_THRESH)
    radar.add_field_like(phidpfield, 'meteoMask', meteoMask)

    # remove phi_dp wrap-around:
    #print 'unwrapping phidp ...'
    (phidp_unwrap) = kdpfun.unwrap_phidp(elev, rays, bins, meteoMask, phidp)

    # remove non-meteo data / filter phi_dp:
    #print 'cleaning / filtering phi_dp ...'
    (phidp_meteo) = kdpfun.clean_phidp(elev, rays, bins, phidp_unwrap, meteoMask,
                                FILTER_BINS_1)

    # generate non-rain mask:
    #print 'generating non-rain mask ...'
    (rainMask) = kdpfun.generate_rain_mask(elev, rays, bins, rhohv, RAIN_THRESH)

    # remove non-rain data / filter phi_dp:
    #print 'removing non-rain components from phidp ...'
    (phidp_rain) = kdpfun.clean_phidp(elev, rays, bins, phidp_meteo, rainMask,
                               FILTER_BINS_2)

    # smooth phi_dp (twice):
    #print 'smoothing phi_dp ...'
    (phidp_smooth) = kdpfun.smooth_data(elev, rays, bins, phidp_rain, SMOOTH_BINS_1)
    (phidp_smooth) = kdpfun.smooth_data(elev, rays, bins, phidp_smooth, SMOOTH_BINS_2)
    radar.add_field_like(phidpfield, 'sPhiDP', phidp_smooth)

    # calculate kdpfun:
    #print 'calculating kdpfun ...'
    (kdp) = kdpfun.calc_kdp_v3(elev, rays, bins, binlength, phidp_smooth, rainMask)

    radar.add_field_like(phidpfield, 'KDP_UKMO', np.ma.masked_array(data=kdp,
                                                             mask=phidp.mask))


def remove_nearest_bins(radar):
	# clean up the nearest beans influenced by sidelobs len(radar.fields[a_field]['data'].shape) > 1
	ranges = radar.range['data']
	for a_field in radar.fields:
		if (a_field != "scan_altitude"):
			radar.fields[a_field]['data'][:,np.where(ranges<400)[0]] = np.nan

def shift_ppi(radar, field_list):
	the_list_of_elevations = np.unique(radar.elevation['data'])

	for elevation in the_list_of_elevations:

		sweep = int(np.where(radar.elevation['data'] == elevation)[0][0] / (int(radar.nrays / radar.nsweeps)))
		sweep_ind = (radar.sweep_start_ray_index['data'][sweep], radar.sweep_end_ray_index['data'][sweep])
		elevation_data = np.where(radar.elevation['data'] == elevation)[0]
		original_azimuthes = radar.azimuth['data'][elevation_data].astype(int)
		north_is_at_position = np.where(original_azimuthes == 0)[0]
		radar.azimuth['data'][elevation_data] = np.roll(original_azimuthes, -north_is_at_position)
		original_time = radar.time['data'][elevation_data]
		radar.time['data'][elevation_data] = np.roll(original_time, -north_is_at_position)

		for field in field_list:
			field_data = radar.fields[field]['data'][sweep_ind[0]:sweep_ind[1]]
			radar.fields[field]['data'][sweep_ind[0]:sweep_ind[1]] = np.roll(field_data, -north_is_at_position, axis=0)

	return radar

def preprocessing(radar, vp_mode):
    """Do preprocessing"""
    try:
        add_beam_height(radar)
    except:
        print ('Beam height failed')
        raise

    if vp_mode == 'QVP': remove_nearest_bins(radar)

    if 'KDP' in radar.fields.keys():
        if 'RhoHV' in radar.fields.keys() and 'uPhiDP' in radar.fields.keys():
            try:
                kdp_ukmo(radar)
            except:
                print ('kdp_ukmo failed')
                raise
        elif 'RhoHV' in radar.fields.keys() and 'uPhiDP' not in radar.fields.keys() and 'PhiDP' in radar.fields.keys():
            try:
                kdp_ukmo(radar,phidpfield='PhiDP')
            except:
                print ('kdp_ukmo failed')
                raise
        else:
            print ('kdp_ukmo failed')

    psidp_field = 0.632 * (copy.deepcopy(radar.fields['ZDR']['data'])) ** 1.71

    radar.add_field_like('PhiDP', 'uPsiDP', psidp_field)

def estimate_noise_level(radar, height_limit=17500, range_limit=50):
    # Mask at the data fill level (MO files)
    dBZ = np.ma.masked_less_equal(radar.fields['reflectivity']['data'],-32)
    
    # Mask high value echoes (unlikely to be noise)
    dBZ_m = np.ma.masked_greater_equal(dBZ,30)
    
    # Define a noise location mask based on altitude and range (mainly height)
    mask = np.logical_and(radar.gate_altitude['data']>height_limit,
                          np.tile(radar.range['data'][:]/1000.0,radar.nrays).reshape(radar.nrays,radar.ngates)>range_limit)
    
    mask = np.logical_and(mask,
                          ~dBZ_m.mask)
        
    # Remove range correction
    dBZ_Plin = dBZ-20*np.log10(radar.range['data']/1000.0)
    # Mask reflectivity with range correction removed
    dBZ_Plin_M = np.where(mask, dBZ_Plin, np.nan)
    noise_estimate = np.round(np.nanmean(np.nanmean(dBZ_Plin_M.data,axis=0)),2)

    lin_noise_estimate = np.power(10,0.1*noise_estimate)
    signal = np.power(10,0.1*dBZ_Plin)-lin_noise_estimate
    
    SNRH = 10*np.log10(signal)-(noise_estimate)
    return(noise_estimate, SNRH)