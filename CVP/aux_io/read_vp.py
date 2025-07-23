from .vp import VerticalProfile
from netCDF4 import Dataset, date2num, num2date
    
def read_vertical_profile(vp_pathname):
    try:    
        vp_data=Dataset(vp_pathname, "r", format="NETCDF4")
        global_attrs={}
        dates=[]
        heights=[]
        max_counts=[]
        longitudes=[]
        latitudes=[]
        print(vp_data)
        for attrname in vp_data.ncattrs():
            this_attr=getattr(vp_data, attrname)
            if attrname=='longitude':
                longitudes=[this_attr]
            elif attrname=='latitude':
                latitudes=[this_attr]
            else:
                global_attrs[attrname]=this_attr
        for var in vp_data.variables:
            if var=='Time':
                time_coord=vp_data.variables[var]
                times=time_coord[:]
                dates=num2date(times,time_coord.units,time_coord.calendar)
            elif var=='Height':
                heights_coord=vp_data.variables[var]
                heights=heights_coord[:]
                height_units=heights_coord.units
            elif var=='Max_Counts':
                max_counts=vp_data.variables[var]
            else:
                print('unexpected variable in VP', var)
        if len(dates)==0 or len(heights)==0:
            print('No times or heights found in VP!')

        long_names={}
        short_names={}
        units={}
        means={}
        stddevs={}
        counts={}
        for field in vp_data.groups:
            var_data=vp_data[field]
            for attrname in var_data.ncattrs():
                this_attr=getattr(var_data, attrname)
                if attrname=='long_name':
                    long_names[field]=this_attr
                elif attrname=='standard_name':
                    short_names[field]=this_attr
                elif attrname=='units':
                    units[field]=this_attr
                else:
                    print('unexpected attribute', attrname, 'in group', field)
            for dat in var_data.variables:
                if dat=='Means':
                    means[field]=var_data.variables[dat][:]
                elif dat=='StdDevs':
                    stddevs[field]=var_data.variables[dat][:]
                elif dat=='Counts':
                    counts[field]=var_data.variables[dat][:]
                else:
                    print('unexpected variable', dat, 'in group', field)
        vp=VerticalProfile(len(dates), heights, vp_data.groups.keys(), long_names, short_names, units, global_attrs, vp_pathname)
        vp.times=dates
        vp.means=means
        vp.stddevs=stddevs
        vp.counts=counts
        vp.lons=longitudes
        vp.lats=latitudes
        
    except OSError as err:
        print(vp_pathname, err)
        vp=None
    return vp
