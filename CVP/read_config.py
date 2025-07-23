#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
read_config function to read the config file passed to vp_extraction

Author: Julia Crook Jun 2025
:copyright: © 2022 University of Leeds.
:license: BSD3

"""

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

