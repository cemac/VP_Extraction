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

#def attenuation_correction(radar):
#	ALPHA = 0.27
#	b = 0.78
#	ZH_AH_DELTAPHI_THRESHOLD = 5
#	AH_FIELD = "Specific_Attenuation_H"
#	Z_IAH_FIELD = "dBZ_ac"
#	ZDR_IAH_FIELD = "ZDR_ac"
#	ZH_AH_FIELD = "Zh_Ah"
#	REFL_METEO = "dBZ"
#	ZDR_METEO = "ZDR"
#	ZH_AH_THR_FIELD = "Zh_Ah_Threshold"
#
#	temperature = copy.deepcopy(radar.fields["temperature_2"]['data'])
#	smooth_phi = copy.deepcopy(radar.fields["sPhiDP"]['data'])
#	meteoMask = copy.deepcopy(radar.fields["meteoMask"]['data'])
#
#	attenuation_mask = np.where(np.logical_or(meteoMask,temperature<273.15),True,False)
#
#	Ah, delta_phi = nrt_attenuation.specific_attenuation_single_segment_ZPHI(radar,
#                                                                            ALPHA,
#                                                                            b,
#                                                                            "sPhiDP",
#                                                                            REFL_METEO,
#                                                                            attenuation_mask,
#                                                                            20, 20, 5, 500,
#                                                                            0, 10, 5, 900)
#
#	delta_phi_array = np.repeat(delta_phi, radar.ngates).reshape(delta_phi.shape[0], radar.ngates)
#
#	radar.add_field_like('PhiDP', 'delta_phi', np.ma.masked_invalid(delta_phi_array))
#	radar.add_field_like('dBuZ', AH_FIELD, np.ma.masked_invalid(Ah))
#
#   # Add attenuation correction to new reflectivity field
#	Z_IAH = nrt_attenuation.correct_Zh_with_Ah(radar.fields[REFL_METEO]['data'],
#                                               radar.fields[AH_FIELD]['data'],
#                                               int(radar.range['meters_between_gates']) / 1000.0)
#
#	radar.add_field_like('dBuZ', Z_IAH_FIELD, Z_IAH)
#
#	ZDR_IAH = nrt_attenuation.correct_zdr_with_Ah(radar.fields[ZDR_METEO]['data'],
#                                                  radar.fields[AH_FIELD]['data'],
#                                                  int(radar.range['meters_between_gates']) / 1000.0,
#                                                  0.14)
#
#	radar.add_field_like('ZDR', ZDR_IAH_FIELD, ZDR_IAH)
#
#	return radar


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

    if vp_mode == 'qvp': remove_nearest_bins(radar)

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
