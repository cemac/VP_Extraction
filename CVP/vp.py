""" 
 Code to build, hold and access vertical profile data
 This only handles static radars where centre lat and lon of the radar should not move and altitudes
 should be the same for each time

 Authors:
 * Written by Julia Crook, September 2024

""" 
import numpy as np
from netCDF4 import Dataset, date2num, num2date
import datetime as dt

#class to hold data in one VP file
class VerticalProfile():
    def __init__(self, ntimes, heights, field_list, long_names, short_names, units, global_attrs, output_file):
        self.times=[dt.datetime(1900,1,1)]*ntimes  # these will be stored as datetimes
        self.heights=heights
        self.output_file=output_file
        self.means={}
        self.stddevs={}
        self.counts={}
        self.long_names={}
        self.standard_names={}
        self.units={}
        self.lons=np.zeros(ntimes)
        self.lats=np.zeros(ntimes)
        self.global_attrs=global_attrs
        nheights=len(heights)
        for field in field_list:
            self.means[field]=np.zeros((nheights, ntimes))*np.nan             
            self.stddevs[field]=np.zeros((nheights, ntimes))*np.nan             
            self.counts[field]=np.zeros((nheights, ntimes), int)
            if field in units.keys():
                self.long_names[field]=long_names[field]
                self.standard_names[field]=short_names[field]
                self.units[field]=units[field]
            else:
                raise ValueError('VerticalProfile creation: units does not contain '+field)
        self.max_counts=np.zeros((nheights, ntimes), int)

    def set_lat_lon_and_bounds(self, centre_lat_lon, lat_lon_bounds):
        self.centre_lat_lon=centre_lat_lon
        self.lat_lon_bounds=lat_lon_bounds
        
    def add_data_for_time(self, field, tix, timeofsweep, mean_values, std_values, counts, lon, lat):
        self.means[field][:,tix]=mean_values           
        self.stddevs[field][:,tix]=std_values     
        self.counts[field][:,tix]=counts
        self.times[tix]=timeofsweep
        self.lons[tix]=lon
        self.lats[tix]=lat
        
        
    def get_field(self, field):
        if field in self.means.keys():
            return (self.means[field], self.stddevs[field], self.counts[field], self.long_names[field], self.standard_names[field], self.units[field])
        else:
            raise ValueError('no such field '+field)
        
    def print(self):
        deltat=self.times[1]-self.times[0]
        print(self.global_attrs['profile_type'], self.output_file, '\ndates:',self.times[0], 'to', self.times[-1], 'every', deltat)
        print(len(self.heights),'heights:',self.heights)
        for key in self.means.keys():
           print('\t'+key+' '+self.standard_names[key])
           print('\t\tCounts min and max:',np.nanmin(self.counts[key], axis=1), np.nanmax(self.counts[key],axis=1))
           if np.nanmax(self.counts[key])>0:
               print('\t\tMeans min and max:',np.nanmin(self.means[key]), np.nanmax(self.means[key]))
               print('\t\tStddevs min and max:',np.nanmin(self.stddevs[key]), np.nanmax(self.stddevs[key]))

    def output_netcdf(self, verbose=False):

        # dataset of the output
        output = Dataset(self.output_file, 'w', format='NETCDF4')

        # Create dimensions
        output.createDimension('Height', self.heights.shape[0])
        output.createDimension('Time', len(self.times))
        # Create coordinate variables for 2 dimensions
        times = output.createVariable('Time', np.float64, ('Time',))
        heights = output.createVariable('Height', np.float32, ('Height',))

        time_units = 'seconds since {}-01-01T00:00:00Z'.format(self.times[0].year)
        number_times = [date2num(sweeptime, time_units, calendar='gregorian') for sweeptime in self.times]

        ## Add values to the dimension variables
        heights[:] = self.heights
        times[:] = np.array(number_times)
        times.units = time_units
        times.calendar = "gregorian"
        heights.units = "metres above sea level"

        # lons and lats are expected to be the same so use a global attribute to hold means of lon and lat
        output.longitude=np.mean(self.lons)
        output.latitude=np.mean(self.lats)
        # are the max_counts in each height always the same in every time?
        ncounts=np.asarray([len(count) for count in np.unique(self.max_counts, axis=1)])
        ix=np.where(ncounts>1)
        if len(ix[0])>0:
            # we need to create a 2d max count
            max_counts = output.createVariable('Max_Counts', np.int16, ('Height', 'Time'))
            max_counts[:]=self.max_counts
        else:
            # we can create a 1d max count as same value for each time
            max_counts = output.createVariable('Max_Counts', np.int16, ('Height',))
            max_counts[:]=self.max_counts[:,0]
            
        # Create the field variables
        # for each variable in the list we create several fields
        for field in self.means.keys():
            group_data_construct_means = '/'+field+'/Means'
            group_data_construct_deviations = '/' + field + '/StdDevs'
            group_data_construct_counts = '/' + field + '/Counts'
            temp_means = output.createVariable(group_data_construct_means, np.float32, ('Height','Time'))
            temp_stds = output.createVariable(group_data_construct_deviations, np.float32, ('Height', 'Time'))
            temp_counts = output.createVariable(group_data_construct_counts, np.int16, ('Height', 'Time'))
            temp_means[:] = self.means[field]
            temp_stds[:] = self.stddevs[field]
            temp_counts[:] = self.counts[field]
            output[field].units = self.units[field]
            output[field].long_name = self.long_names[field]
            output[field].standard_name = self.standard_names[field]

        for key in self.global_attrs.keys():
            # add global attribute
            setattr(output, key, self.global_attrs[key])
            
        output.close()

        if verbose:
            print('New netcdf file created at {}'.format(self.output_file))
