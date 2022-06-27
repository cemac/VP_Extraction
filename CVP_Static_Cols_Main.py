"""
Routines for Writing CVP Cols from UKMO NIMROD Files
"""

cvp_script="/gws/smf/j04/ncas_radar/rrniii/BioDAR/CVP_Extraction/cvp_extraction_static.py"
#log output location
output="/gws/smf/j04/ncas_radar/rrniii/BioDAR/CVP_Extraction/Output/"
    
from datetime import datetime, timedelta
import os
import sys

#function for getting list of dates
def date_range(start, end):
    delta = end - start  # as timedelta
    start_list = [start + timedelta(days=i) for i in range(delta.days + 1)]
    end_list   = [start + timedelta(days=j, hours=23, minutes=59, seconds=59) for j in range(delta.days + 1)]
    return start_list, end_list
        
### Main Program ###
#Define Sites with Nested Dictionary
cvp_sites = {
         0: {'col_pos_name': "rothamsted",            'col_lat' : 51.806908, 'col_long' : -0.3609735, 'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         1: {'col_pos_name': "rothamsted_fixed_15km", 'col_lat' : 51.78928,  'col_long' : -0.38672,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         2: {'col_pos_name': "rothamsted_fixed_25km", 'col_lat' : 51.85451,  'col_long' : -0.28863,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         3: {'col_pos_name': "rothamsted_fixed_35km", 'col_lat' : 51.92049,  'col_long' : -0.18684,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         4: {'col_pos_name': "silwood",               'col_lat' : 51.409459, 'col_long' : -0.6426696, 'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         5: {'col_pos_name': "silwood_fixed_15km",    'col_lat' : 51.55378,  'col_long' : -0.57675,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         6: {'col_pos_name': "silwood_fixed_25km",    'col_lat' : 51.46861,  'col_long' : -0.60742,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         7: {'col_pos_name': "silwood_fixed_35km",    'col_lat' : 51.37699,  'col_long' : -0.65068,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
         8: {'col_pos_name': "wye",                   'col_lat' : 51.185132, 'col_long' : 0.9448683,  'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
         9: {'col_pos_name': "wye_fixed_15km",        'col_lat' : 51.22918,  'col_long' : 0.79166,    'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
        10: {'col_pos_name': "wye_fixed_25km",        'col_lat' : 51.18637,  'col_long' : 0.9314,     'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
        11: {'col_pos_name': "wye_fixed_35km",        'col_lat' : 51.14978,  'col_long' : 1.05019,    'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
        12: {'col_pos_name': "hereford",              'col_lat' : 52.124081, 'col_long' : -2.6382625, 'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
        13: {'col_pos_name': "hereford_fixed_15km",   'col_lat' : 52.25514,  'col_long' : -2.6118,    'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
        14: {'col_pos_name': "hereford_fixed_25km",   'col_lat' : 52.17268,  'col_long' : -2.62416,   'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
        15: {'col_pos_name': "hereford_fixed_35km",   'col_lat' : 52.08205,  'col_long' : -2.64367,   'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
        16: {'col_pos_name': "preston",               'col_lat' : 53.854688, 'col_long' : -2.7646454, 'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
        17: {'col_pos_name': "preston_fixed_15km",    'col_lat' : 53.80162,  'col_long' : -2.51104,   'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
        18: {'col_pos_name': "preston_fixed_25km",    'col_lat' : 53.82787,  'col_long' : -2.65404,   'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
        19: {'col_pos_name': "preston_fixed_35km",    'col_lat' : 53.86671,  'col_long' : -2.79909,   'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
        20: {'col_pos_name': "starcross",             'col_lat' : 50.629506, 'col_long' : -3.4548264, 'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
        21: {'col_pos_name': "starcross_fixed_15km",  'col_lat' : 50.82359,  'col_long' : -3.45459,   'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
        22: {'col_pos_name': "starcross_fixed_25km",  'col_lat' : 50.73279,  'col_long' : -3.45605,   'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
        23: {'col_pos_name': "starcross_fixed_55km",  'col_lat' : 50.63798,  'col_long' : -3.45654,   'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
        } 

#Define Start and End
start_date = datetime(2019, 1, 2)
end_date = datetime(2019, 12, 31)
    
[start_list, end_list] = date_range(start_date, end_date)
    
print("Overall start date is: ", start_list[0].strftime("%Y%m%d"+"T"+"%H%M%S"))
print("Overall end date is:   ", end_list[-1].strftime("%Y%m%d"+"T"+"%H%M%S"))

index = sys.argv[1]

d = index//len(cvp_sites)
s = index%len(cvp_sites)
        
year=str(start_list[d].year)
start = start_list[d].strftime("%Y%m%d"+"T"+"%H%M%S")
end = end_list[d].strftime("%Y%m%d"+"T"+"%H%M%S")

col_pos_name = cvp_sites[s]['col_pos_name']
col_long = cvp_sites[s]['col_long']
col_lat = cvp_sites[s]['col_lat']
radar = cvp_sites[s]['radar']

for r in range(0, len(cvp_sites[s]['col_radius'])):

    col_radius = cvp_sites[s]['col_radius'][r]

    cvp_input_data_dir="/gws/smf/j07/ncas_radar/data/ukmo-nimrod/raw_h5_data/single-site/"+radar+"/"+year+"/"
    print(cvp_input_data_dir)

    cvp_out_data_dir="/gws/smf/j07/ncas_radar/data/ukmo-nimrod/raw_cvp_data/"+radar+"/"+col_pos_name+"/"+str(col_radius).replace(".","_")+"/"+year+"/"
    print(cvp_out_data_dir)

    file_check=cvp_out_data_dir+col_pos_name+"_"+f'{col_radius:.1f}'+"km_"+str(start_list[d].year)+f'{start_list[d].month:02}'+f'{start_list[d].day:02}'+".nc"
    print(file_check)
    if os.path.exists(file_check):
        print("Done Already: "+file_check)
        with open(output+'CVP_Cols_Done.txt', 'a') as log:
            log.write("Done Already: "+file_check+ '\n')
    else:
        job_name= radar+"_"+start+"_"+col_pos_name+"_"+str(col_radius).replace(".","_")
        print(job_name)
        with open(output+'CVP_Cols_To_do.txt', 'a') as log:
            log.write("Sent to Lotus: "+file_check+ '\n')
        python_run_cmd='python '+str(cvp_script)+' '+str(cvp_input_data_dir)+' '+str(cvp_out_data_dir)+' -d -v -m -s '+str(start)+' -e '+str(end)+' -r '+str(col_radius)+' -a '+str(col_lat)+' -o '+str(col_long)+' -p '+str(col_pos_name)+'"'
                    
        #print python command for debugging
        #print('python '+ str(cvp_script)+' '+str(cvp_input_data_dir)+' '+str(cvp_out_data_dir) +' -d -v -m -s ' +str(start)+ ' -e ' +str(end)+' -r ' + str(col_radius) +' -a ' +str(col_lat) + ' -o ' + str(col_long) +' -p ' + str(col_pos_name))
                    
        os.system(python_run_cmd)
                   
                    
#queues --account=short4hr --partition=short-serial
#--partition=tes
