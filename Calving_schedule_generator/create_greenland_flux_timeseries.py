
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_mankoff_discharge_data(data_dir):
    ds = nc4.Dataset(os.path.join(data_dir,'gate.nc'))
    discharge = ds.variables['discharge'][:,:]
    time = ds.variables['time'][:]
    gates = ds.variables['gate'][:]
    name_Mouginot = list(ds.variables['name_Mouginot'][:])
    ds.close()

    dec_yrs = np.zeros((len(time),))
    for t in range(len(time)):
        date = datetime.datetime(1986,4,14) + datetime.timedelta(days=int(time[t]))
        dec_yrs[t] = YMD_to_DecYr(date.year, date.month, date.day)

    all_timeseries = []

    # for gate in range(len(name_Mouginot)):
    #     if name_Mouginot[gate] in glacier_names:
    #         print(gates[gate])
    #         timeseries = np.column_stack([time,discharge[gate,:]])
    #         plt.plot(timeseries[:,0], timeseries[:,1])
    #         plt.show()
    #         all_timeseries.append(timeseries)

    for gate in gates:
        index = np.where(gates==gate)[0][0]
        print('Reading '+name_Mouginot[index])
        # print(index,np.shape(dec_yrs))
        timeseries = np.column_stack([dec_yrs,discharge[index,:]])
        # plt.plot(timeseries[:,0], timeseries[:,1])
        # plt.show()
        all_timeseries.append(timeseries)
        name_Mouginot[index] = str(gate) + '_' + str(name_Mouginot[index])

    return name_Mouginot, all_timeseries

def save_discharge_timeseries_as_nc(project_dir,glacier_names,all_timeseries):

    if not os.path.exists(os.path.join(project_dir,'Data','Glacier')):
        os.makedirs(os.path.join(project_dir,'Data','Glacier'))
    if os.path.exists(os.path.join(project_dir,'Data','Glacier','Ice Flux Timeseries.nc')):
        os.remove(os.path.join(project_dir,'Data','Glacier','Ice Flux Timeseries.nc'))

    output_file = os.path.join(project_dir,'Data','Glacier','Ice Flux Timeseries.nc')
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('dec_yr',np.shape(all_timeseries[0])[0])

    d = ds.createVariable('dec_yr','f4',('dec_yr',))
    d[:] = all_timeseries[0][:,0]

    for g in range(len(glacier_names)):
        v = ds.createVariable(glacier_names[g],'f4',('dec_yr',))
        v[:] = all_timeseries[g][:,1]

    ds.close()


def camelize(x: str): # To transform all capitalized to camel case
    parts=x.split('_')
    parts = [part.capitalize() for part in parts]
    return '_'.join(parts)

project_dir = '/Volumes/T7/IceBerg'

data_dir = '/Volumes/T7/IceBerg/flux_data'

# mouginot_glacier_names = ['UPERNAVIK_ISSTROM_SS',
#                  'UPERNAVIK_ISSTROM_S',
#                  'UPERNAVIK_ISSTROM_C',
#                  'UPERNAVIK_ISSTROM_N']
#
# gate_numbers = [129,130,131,132,133]
#
# glacier_names = ['Upernavik_Isstrom_NW',
#                  'Upernavik_Isstrom_N',
#                  'Upernavik_Isstrom_C',
#                  'Upernavik_Isstrom_S',
#                  'Upernavik_Isstrom_SS']

glacier_names, all_timeseries = read_mankoff_discharge_data(data_dir)

# glacier_names = [camelize(name) for name in glacier_names]

save_discharge_timeseries_as_nc(project_dir,glacier_names,all_timeseries)

