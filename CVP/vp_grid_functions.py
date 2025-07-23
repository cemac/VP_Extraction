# CVP grid related functions for constructing the grid that don't need radar data
# Author: Julia Crook, June 2025

import numpy as np
import math

#-----------------------------------------------------------
# find the new lat and lon in degrees from lat1, lon1 (in degrees)
# that is distance km away along the azimuth
#-----------------------------------------------------------
def get_new_lat_lon(lat1, lon1, distance, azimuth):
    # from https://www.movable-type.co.uk/scripts/latlong.html
    Re=6378 # radius of the earth in km
    phi1=np.radians(lat1)
    lambda1=np.radians(lon1)
    azimuth_rad=np.radians(azimuth)
    phi2 = math.asin( math.sin(phi1)*math.cos(distance/Re) + math.cos(phi1)*math.sin(distance/Re)*math.cos(azimuth_rad) )
    lambda2 = lambda1 + math.atan2( math.sin(azimuth_rad)*math.sin(distance/Re)*math.cos(phi1), 
                                    math.cos(distance/Re) - math.sin(phi1) * math.sin(phi2) )
    return np.degrees(phi2), np.degrees(lambda2)

#----------------------------------------------------------------------------
# Calculate the lats and lons of the centre of each of the nrows x ncolumns grids
# First find azimuths and horizontal radii for the centre of each 
# We only go out to max_radius and each cvp has a radius of col_radius
# This allows us to calculate nrows and ncolumns and the positions of the centres of each cvp
# Also calculate the lat lon boundaries of the CVPs for square CVPs
# Inputs:
#    centre_lat, centre_lon - the centre of teh radar from which we construct the grid
#    col_radius - the radius in km of each CVP
#    max_radius - how far to go out from the centre
# ----------------------------------------------------------------------------
def get_cvp_grid_lat_lons(centre_lat, centre_lon, col_radius, max_radius=30):
    ncolumns=nrows=int(max_radius/col_radius)
    ngrids=ncolumns*nrows
    # distance between cvp centres is 2*col_radius left most cvps start at max_radius km distance west from centre
    Xs=np.arange(ncolumns)*2*col_radius-max_radius+col_radius
    radii=np.zeros(ngrids)
    azimuths=np.zeros(ngrids)
    gid=0
    for r in range(nrows):
        y=Xs[r]
        thetas=np.asarray([math.atan2(y,Xs[i]) for i in range(ncolumns)])
        radii[gid:gid+ncolumns]=np.asarray([Xs[i]/math.cos(thetas[i]) for i in range(ncolumns)])
        # theta is measured anticlockwise from east, azimuth is measured clockise from north and is in degrees
        azimuths[gid:gid+ncolumns]=(np.pi/2-thetas)*180/np.pi
        gid+=ncolumns        
    # azimuths are always positive in radar
    ix=np.where(azimuths<0)
    azimuths[ix[0]]+=360

    # now convert these to lats and lons
    grid_lat_lons=np.asarray([get_new_lat_lon(centre_lat, centre_lon, radii[gid], azimuths[gid]) for gid in range(ngrids)])
    # to find the boundaries just use half way between the centres in terms of lats and lons
    # we could use get_specific_cvp_lat_lon_bounds but this tends to overlap the boundaries slightly
    grid_lat_lons_reshape=grid_lat_lons.reshape(nrows, ncolumns, 2)
    grid_bounds=np.zeros((nrows, ncolumns, 4, 2))
    right_deltas=np.zeros((nrows, ncolumns, 2))
    left_deltas=np.zeros((nrows, ncolumns, 2))
    above_deltas=np.zeros((nrows, ncolumns, 2))
    below_deltas=np.zeros((nrows, ncolumns, 2))
    for c in range(ncolumns):
        if c==ncolumns-1:
            # use half the difference to the left
            right_deltas[:,c,:]=(grid_lat_lons_reshape[:,c,:]-grid_lat_lons_reshape[:,c-1,:])/2
        else:
            # use half the difference to the right
            right_deltas[:,c,:]=(grid_lat_lons_reshape[:,c+1,:]-grid_lat_lons_reshape[:,c,:])/2
        if c==0:
            # use right_deltas
            left_deltas[:,c,:]=right_deltas[:,c,:]
        else:
            # use half the difference to the left
            left_deltas[:,c,:]=(grid_lat_lons_reshape[:,c,:]-grid_lat_lons_reshape[:,c-1,:])/2

    for r in range(nrows):
        if r==nrows-1:
            # use half the difference to below row
            above_deltas[r,:,:]=(grid_lat_lons_reshape[r,:,:]-grid_lat_lons_reshape[r-1,:,:])/2
        else:
            # use half the difference to above row
            above_deltas[r,:,:]=(grid_lat_lons_reshape[r+1,:,:]-grid_lat_lons_reshape[r,:,:])/2
        if r==0:
            # use the above_deltas
            below_deltas[r,:,:]=above_deltas[r,:,:]
        else:
            # use half the difference to below row
            below_deltas[r,:,:]=(grid_lat_lons_reshape[r,:,:]-grid_lat_lons_reshape[r-1,:,:])/2
    # top left
    grid_bounds[:,:,0,:]=grid_lat_lons_reshape-left_deltas+above_deltas
    # top right
    grid_bounds[:,:,1,:]=grid_lat_lons_reshape+right_deltas+above_deltas
    # bottom right
    grid_bounds[:,:,2,:]=grid_lat_lons_reshape+right_deltas-below_deltas
    # bottom left
    grid_bounds[:,:,3,:]=grid_lat_lons_reshape-left_deltas-below_deltas
    grid_bounds=grid_bounds.reshape(ngrids,4,2)

    return grid_lat_lons, grid_bounds
    
def get_specific_cvp_lat_lon_bounds(centre_lat, centre_lon, col_radius):
    hypotenuse=col_radius/math.cos(np.radians(45))
    lat_lon_top_left=get_new_lat_lon(centre_lat, centre_lon, hypotenuse, 315)
    lat_lon_top_right=get_new_lat_lon(centre_lat, centre_lon, hypotenuse, 45)
    lat_lon_bottom_left=get_new_lat_lon(centre_lat, centre_lon, hypotenuse, 225)
    lat_lon_bottom_right=get_new_lat_lon(centre_lat, centre_lon, hypotenuse, 135)

    return np.asarray([lat_lon_top_left, lat_lon_top_right, lat_lon_bottom_right,lat_lon_bottom_left])
    
