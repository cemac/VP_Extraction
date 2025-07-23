#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Vertical profile extraction

Extract QVP/CVP from radar scan data and store the profile in a
specified format.

This script was developed by CEMAC as part of the Drivers and
Repercussions of UK Insect Declines (DRUID) project (NE/V006916/1).

Authors:
 * Written by David Dufton, Nov. 2016
 * Adapted by M. Lukach, May 2018
 * Adapted by R. Neely, March 2022
 * Adapted by T.D. James <t.d.james1@leeds.ac.uk>, June 2022
 * Adapted by J.A. Crook June 2025 to read each file and do preprocessing and then loop round each cvp or qvp 

:copyright: © 2022 University of Leeds.
:license: BSD3

"""
from __future__ import (absolute_import, division, print_function)
#from six.moves import (filter, input, map, range, zip)

import sys
import os
import glob
import datetime as dt
import argparse
import itertools
import numpy as np
from dateutil.parser import parse as dateparse
import h5py as h5

import vp_functions
import vp_io
from vp import *
import read_config
import vp_grid_functions
import pdb

def parse_args():
    formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter)

    parser.add_argument("-r", "--radar",
                        dest='radar_name',
                        help='''radar name''')
    parser.add_argument("-c", "--cfg_file",
                        dest='config_file',
                        help='''config file defining how to extract (qvp/cvp) and in/out directories''')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help='''Print messages about the program execution
                        to the console (stdout)''')

    parser.add_argument("-t", "--time",
                        dest="timestamp",
                        help='''timestamp of raw radar file:
                        "YYYYMMDD" as used in raw radar directories/filenames''')

    
    args = parser.parse_args()
    return args


def get_input_folder_glob_spec(input_dir,
                               profile_type,
                               vertical,
                               timestamp,
                               met_office,
                               verbose=False):
    '''
    Get pattern for matching files in input folder.

    The matching pattern depends on profile type (vertical scan for QVP) and source of radar scan data.
    '''

    # Path to the directory with the radar scans
    folder_with_files = os.path.normpath(input_dir)
    if not os.path.isdir(folder_with_files):
        folder_with_files = os.path.dirname(folder_with_files)
    if verbose:
        print("Path to input directory is", folder_with_files)

    # Handle different options for input directory structure
    if met_office:
        folder_glob_spec = '{}/{}_*.h5'.format(folder_with_files, timestamp)
    else:
        folder_glob_spec = '{}/*{}*.nc'.format(folder_with_files, timestamp)
    if verbose:
        print(("Input folder glob spec is {}").format(folder_glob_spec))

    return folder_glob_spec


def get_file_list(input_dir,
                  profile_type,
                  vertical,
                  timestamp,
                  met_office,
                  verbose=False):
    '''
    Get list of files for profile extraction.

    '''

    # Get file matching pattern
    folder_glob_spec = get_input_folder_glob_spec(input_dir,
                                                  profile_type,
                                                  vertical,
                                                  timestamp,
                                                  met_office,
                                                  verbose)

    # select all files available in the directory
    file_list = glob.glob(folder_glob_spec)
    file_list.sort()

    if verbose:
        if file_list:
            print('First file is:', file_list[0])
            print('Last file is:', file_list[-1])
            print("There are ", len(file_list), " files")

    return file_list

def main():

    #### Input management

    args = parse_args()
    config_file=args.config_file
    config = read_config.read_config(args.config_file)
    # get the timestamp as a datetime
    timestamp=args.timestamp
    t_datetime = dateparse(timestamp)
    
    verbose = args.verbose

    # True if the data are from the met office radar
    met_office = config['MET_OFFICE']

    # Profile type (QVP or CVP)
    profile_type = config['PROFILE_TYPE']
    
    input_dir=t_datetime.strftime(config['DATA_INPUT']) # add defined part of the timestamp (could be just year or whole timestamp)
    # Check if input directory exists
    if not os.path.exists(input_dir):
        err_msg = "Input dir {0} does not exist\n"
        err_msg = err_msg.format(input_dir)
        raise ValueError(err_msg)
    output_dir=config['DATA_OUTPUT']
    
    if profile_type=='CVP':
        vertical=False # only used for QVPs
        nheights=int(config['MAX_H']/config['H_STEP'])
        equidistant_alt = np.linspace((config['MIN_H'] + config['H_STEP']/2), (config['MAX_H'] - config['H_STEP']/2), num=nheights)
        equidistant_bound = np.linspace((config['MIN_H']), (config['MAX_H']), num=nheights+1)
        output_dir=output_dir+'{}km/'.format(config['COL_RADIUS'])
        ngrids=0
    else:
        azimuth_exclude=config['AZIMUTHS_TO_EXCLUDE']
        elevations=config['ELEVATIONS']
        vertical=elevations[0]==90 # if we are doing vertical radar sweep this should be the only elevation as the files are else where

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # File list
    file_list = get_file_list(input_dir,
                              profile_type,
                              vertical,
                              args.timestamp,
                              met_office,
                              verbose)

    # read the files one at a time and one time from them at a time
    defaulttime = [dt.datetime(1970,1,1,0,0,0)]
    if(met_office):
        # there should only be one file so get the times from it
        testfile = h5.File(file_list[0], 'r')
        times = list(testfile['lp'].keys())
        file_list = list(itertools.chain.from_iterable(itertools.repeat(x, len(times)) for x in file_list))
        testfile.close()
    else:
        # there will be files with each time for the day
        times = [defaulttime]

    vps=[]
    unit_dict = []
    long_names = []
    short_names = []
    ntimes=len(file_list)
    for f, file_ in enumerate(file_list):
        if (verbose):
            print ("file f {} is {}".format(f,file_))

        # Edited function call so that different call not needed for metoffice=True
        (radar, unit_dict, long_names, short_names) = \
            vp_io.read_file(f, times[f%len(times)], file_list[f], config['LOG_OUTPUT'], config['FIELD_LIST'], unit_dict,
                      long_names, short_names, met_office=met_office, verbose=verbose)

        radar_lat, radar_lon=vp_functions.get_centre_lat_lon_for_radar(radar)
        timeofsweep = num2date(np.nanmean(radar.time['data'][:]),
                                          radar.time['units'],
                                          radar.time['calendar'])
        if verbose:
            print('radar at lat lon', radar_lat, radar_lon)
        nazimuths=int(radar.nrays/radar.nsweeps)
        if profile_type == "QVP":
            for i in range(len(elevations)):
                if f==0:
                    output_file=output_dir+'{}_QVP_{:.1f}deg.nc'.format(args.timestamp, elevations[i])
                    # create a vertical profile object
                    radar_elevations = radar.elevation['data'].reshape((radar.nsweeps,nazimuths))
                    radar_elevations=radar_elevations[:,0] # elevations are same for all azimuths
                    eix=np.argmin(radar_elevations-elevations[i])
                    # getting altitudes here assumes altitudes are always the same
                    if elevations[i] == 90.0:
                        altitudes = radar.range['data']
                    else:
                        altitudes = radar.fields['scan_altitude']['data'].reshape(radar.nsweeps,nazimuths,radar.ngates)
                        altitudes=altitudes[eix,0,:]
                        
                    global_attrs={}
                    global_attrs['radar_name']=args.radar_name
                    global_attrs['profile_type']=profile_type
                    global_attrs['elevation']=elevations[i]
                    global_attrs['from_input']=input_dir
                    azimuth_exc_strs=[str(from_to[0])+'-'+str(from_to[1]) for from_to in azimuth_exclude]
                    azimuth_exc_str=', '.join(azimuth_exc_strs)
                    global_attrs['azimuths_excluded']=azimuth_exc_str
                    vps.append(VerticalProfile(ntimes, altitudes, config['FIELD_LIST'], long_names, short_names, unit_dict, global_attrs, output_file))
                    if verbose:
                        print("Output will be placed in {} for elevation".format(output_file), elevations[i])
                else:
                    if elevations[i] == 90.0:
                        altitudes = radar.range['data']
                    else:
                        altitudes = radar.fields['scan_altitude']['data'].reshape(radar.nsweeps,nazimuths,radar.ngates)
                        altitudes=altitudes[eix,0,:]
                    max_diff=np.amax(abs(vps[i].heights-altitudes))
                    if max_diff>0.1:
                        raise Exception('QVP has different altitudes at different times')
                this_vp=vps[i]
                # Extract QVP
                if verbose:
                    print('extracting QVP for elevation', elevations[i], 'for time', timeofsweep)
                vp_functions.time_height_qvp(radar,
                                             this_vp,
                                             config['FIELD_LIST'],
                                             f,
                                             elevations[i],
                                             azimuth_exclude=azimuth_exclude,
                                             verbose=verbose)
                 

        elif profile_type == 'CVP':
            if config['CREATE_CVP_GRID']:
                if f==0:
                    # get the centre lat lon of the radar and then calculate the lat lons of the grids
                    grid_lat_lons, grid_bounds=vp_grid_functions.get_cvp_grid_lat_lons(radar_lat, radar_lon,
                                                                          config['COL_RADIUS'],config['MAX_RADIUS'])
                    ngrids=grid_lat_lons.shape[0]

                for gid in range(ngrids):
                    if f==0:
                        site_name='{}_GridID_{:03d}'.format(args.radar_name, gid+1)
                        this_output_dir=output_dir+site_name+'/{}/'.format(t_datetime.year)
                        if not os.path.exists(this_output_dir):
                            os.makedirs(this_output_dir)
                        output_file=this_output_dir+'{}_{}km_{}.nc'.format(site_name,
                                                                          config['COL_RADIUS'],
                                                                          args.timestamp)
                        global_attrs={}
                        global_attrs['radar_name']=args.radar_name
                        global_attrs['profile_type']=profile_type
                        global_attrs['from_input']=input_dir
                        vps.append(VerticalProfile(ntimes, equidistant_alt, config['FIELD_LIST'], long_names, short_names, unit_dict, global_attrs, output_file))
                        vps[-1].set_lat_lon_and_bounds(grid_lat_lons[gid,:], grid_bounds[gid,:,:])

                    this_vp=vps[gid]
                    # Extract CVP
                    if verbose:
                        print('extracting CVP grid', gid+1, 'for time', timeofsweep)
                    vp_functions.time_height_cvp(radar,
                                                 this_vp,
                                                 config['COL_RADIUS'],
                                                 config['FIELD_LIST'],
                                                 f,
                                                 equidistant_bound,
                                                 verbose=verbose)
            for item, specific_cvp in enumerate(config['SPECIFIC_CVP']):
                if f==0:
                    site_name='{}_{}'.format(args.radar_name, specific_cvp[2])
                    this_output_dir=output_dir+site_name+'/{}/'.format(t_datetime.year)
                    if not os.path.exists(this_output_dir):
                        os.makedirs(this_output_dir)
                    output_file=this_output_dir+'/{}_{}km_{}.nc'.format(site_name,
                                                                       config['COL_RADIUS'],
                                                                       args.timestamp)
                    global_attrs={}
                    global_attrs['radar_name']=args.radar_name
                    global_attrs['profile_type']=profile_type
                    global_attrs['from_input']=input_dir
                    vps.append(VerticalProfile(ntimes, equidistant_alt, config['FIELD_LIST'], long_names, short_names, unit_dict, global_attrs, output_file))
                    lat_lon_bounds=vp_grid_functions.get_specific_cvp_lat_lon_bounds(specific_cvp[0], specific_cvp[1], config['COL_RADIUS'])
                    vps[-1].set_lat_lon_and_bounds(specific_cvp[:2], lat_lon_bounds)
                this_vp=vps[ngrids+item]
                # Extract CVP
                if verbose:
                    print('extracting specific CVP', specific_cvp[2], 'for time', timeofsweep)
                vp_functions.time_height_cvp(radar,
                                             this_vp,
                                             config['COL_RADIUS'],
                                             config['FIELD_LIST'],
                                             f,
                                             equidistant_bound,
                                             verbose=verbose)
        else:
            raise ValueError("Unrecognised profile type")

    #### Output management

    # Output all vps to NetCDF
    for vp in vps:
        if verbose:
            # lons and lats are expected to be the same 
            print(vp.output_file, 'cvp, lon:', np.mean(vp.lons),'lat:', np.mean(vp.lats), 'complete')

        vp.output_netcdf(verbose=verbose)
    with open(config['LOG_OUTPUT'], 'a') as log:
        log.write(dt.datetime.today().strftime('%Y-%m-%d %H:%M: ')+profile_type+' created for '+input_dir+ '\n')

    # end main()


if __name__ == "__main__":
    main()
