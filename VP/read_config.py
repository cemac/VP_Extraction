#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
read_config function to read the config file passed to vp_extraction

Author: Julia Crook Jun 2025
:copyright: © 2022 University of Leeds.
:license: BSD3

"""

import configparser as ConfParse

def read_config(cfg_file):
    config={}
    specific_cvps=[]
    azimuths_to_exclude=[]
    required_keys=['LOG_OUTPUT', 'DATA_INPUT', 'DATA_OUTPUT', 'MET_OFFICE', 'PROFILE_TYPE', 'FIELD_LIST']
    qvp_keys=['ELEVATIONS', 'AZIMUTHS_TO_EXCLUDE']
    cvp_keys=['COL_RADIUS', 'MAX_RADIUS', 'MIN_H', 'MAX_H', 'H_STEP', 'CREATE_CVP_GRID', 'SPECIFIC_CVP']
    with open(cfg_file,"r") as f:
        lines=f.readlines()
        for line in lines:
            line=line.rstrip()
            wsplit=line.split(' ')
            if wsplit[0] not in required_keys and wsplit[0] not in qvp_keys and wsplit[0] not in cvp_keys and wsplit[0]!='#':
                print('Warning: unknown key in config file', wsplit[0])
            if len(wsplit)<2 or wsplit[0]=='#':
                continue
                
            key=wsplit[0]   
            if key=='SPECIFIC_CVP':
                # each specific cvp has lat lon name
                if len(wsplit)>3:
                    this_cvp=[float(wsplit[1]), float(wsplit[2]), wsplit[3]]
                    specific_cvps.append(this_cvp)
                else:
                    print('no specific cvps defined')
            elif key=='AZIMUTHS_TO_EXCLUDE':
                # each azimuths to exclude has from_azimuth, to_azimuth
                if len(wsplit)>2:
                    azimuths_to_exclude.append([float(wsplit[1]), float(wsplit[2])])
                else:
                    print('no azimuths to exclude defined')
            elif key in ['MIN_H','MAX_H','H_STEP']: # these are integers
                config[wsplit[0]]=int(wsplit[1])
            elif key=='COL_RADIUS': # this is a float
                config[wsplit[0]]=float(wsplit[1])
            elif key=='MAX_RADIUS': # this is a float
                config[wsplit[0]]=float(wsplit[1])
            elif key=='ELEVATIONS':
                config[wsplit[0]]=[float(key) for key in wsplit[1:]]
            else:
                if len(wsplit)>2:
                    config[wsplit[0]]=wsplit[1:]
                else:
                    config[wsplit[0]]=wsplit[1]
                # handle the booleans
                if wsplit[1]=='True':
                    config[wsplit[0]]=True
                elif wsplit[1]=='False':
                    config[wsplit[0]]=False

    config['SPECIFIC_CVP']=specific_cvps
    config['AZIMUTHS_TO_EXCLUDE']=azimuths_to_exclude

    for key in required_keys:
        if key not in config.keys():
            print('Warning: missing required config', key)
            
    if config['PROFILE_TYPE']=='QVP':
        for key in qvp_keys:
            if key not in config.keys():
                print('Warning: missing required QVP config', key)
    elif config['PROFILE_TYPE']=='CVP':
        for key in cvp_keys:
            if key not in config.keys():
                print('Warning: missing required CVP config', key)
    else:
        print('Warning: invalid PROFILE_TYPE in config')
    return config


def read_cvp_config(cfg_file):

    config={}
    required_keys=['LOG_OUTPUT', 'DATA_INPUT', 'DATA_OUTPUT', 'FIELD_LIST']
    cvp_keys=['COL_RADIUS', 'MAX_RADIUS', 'MIN_H', 'MAX_H', 'H_STEP', 'CREATE_CVP_GRID', 'SPECIFIC_CVP']
    
    config_f = ConfParse.ConfigParser(interpolation=None)
    config_f.optionxform = lambda option: option
    config_f.read(cfg_file)

    for key in config_f['FILE_IO']:
        config.update({key: config_f['FILE_IO'][key]})
    
    for key in config_f['CVP_SETTINGS']:
        config.update({key: config_f['CVP_SETTINGS'][key]})

    config['FIELD_LIST'] = []
    for key_ in config_f['FIELD_LIST']:
        config['FIELD_LIST'].append(config_f['FIELD_LIST'][key_])

    config.update({'params_dict': {}})
    for key in config_f['CVP_parameters']:
        config['params_dict'].update({key: config_f.getfloat('CVP_parameters', key)})

    config['CREATE_CVP_GRID'] = config_f.getboolean('CVP_LOCATIONS','CREATE_CVP_GRID')

    for key in ['MIN_H','MAX_H','H_STEP']:
        config[key] = config_f.getint('CVP_SETTINGS',key)

    for key in ['COL_RADIUS','MAX_RADIUS']:
        config[key] = config_f.getfloat('CVP_SETTINGS',key)

    config['SPECIFIC_CVP'] = []
    for key_ in config_f['SPECIFIC_CVPS']:
        loc_string = config_f['SPECIFIC_CVPS'][key_]
        locs = loc_string.split(' ')
        config['SPECIFIC_CVP'].append((float(locs[0]), float(locs[1]), locs[2]))

    for key in required_keys:
        if key not in config.keys():
            print('Warning: missing required config', key)
            
    if config['PROFILE_TYPE']=='CVP':
        for key in cvp_keys:
            if key not in config.keys():
                print('Warning: missing required CVP config', key)
    else:
        print('Warning: invalid PROFILE_TYPE in config')
    return config