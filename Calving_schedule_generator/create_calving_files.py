
import os
import shutil

import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from datetime import datetime, timedelta
import shapefile

def read_grid_geometry(grid_path):

    ds = nc4.Dataset(grid_path)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    return(XC, YC, Depth)

# def read_model_front_locations_from_shp(glacier_names, shapefile_name):
#
#     r = shapefile.Reader(shapefile_name)
#     records = r.records()
#
#     glacier_dict = {}
#     for glacier in glacier_names:
#         rows = []
#         cols = []
#         depths = []
#         for i in range(len(records)):
#             record = list(records[i])
#             if glacier in record:
#                 rows.append(record[0])
#                 cols.append(record[1])
#                 depths.append(int(record[2]))
#         rows = np.array(rows)
#         cols = np.array(cols)
#         depths = np.array(depths)
#         glacier_dict[glacier] = [rows, cols, depths]
#
#     return(glacier_dict)


def mapping_to_row_and_col(coordinate_list, grid_path, n_rows, n_cols):
    model_grid = np.fromfile(os.path.join(grid_path), '>f8').reshape((16, n_rows + 1, n_cols + 1))
    # recreate the grids that will be used in the model
    XC = model_grid[0, :-1, :-1]
    YC = model_grid[1, :-1, :-1]
    coordinate = np.zeros((n_rows, n_cols, 2), dtype=np.float32)
    for row in range(n_rows):
        for col in range(n_cols):
            coordinate[row, col] = (XC[row, col], YC[row, col])

    glacier_map = []
    for loc in coordinate_list:
        glacier_map.append(look_up(loc[0], loc[1], coordinate))
    return glacier_map


def get_distance(x1, y1, x2, y2):
    if x1 != 0 and x2 != 0 and y1 != 0 and y2 != 0:
        return np.hypot(x2 - x1, y2 - y1)
    return 10000000


def look_up(x, y, coordinate):
    min_loc = (-1, -1, 100000)
    for row in range(coordinate.shape[0]):
        for col in range(coordinate.shape[1]):
            distance = get_distance(coordinate[row, col, 0], coordinate[row, col, 1], x, y)
            if distance < min_loc[2]:
                min_loc = (row, col, distance)
    return min_loc[0], min_loc[1]



# Read loc file created by gate_location_extractor.py script
def read_glacier_loc_from_loc_file(loc_file_path: str):
    return np.load(loc_file_path, allow_pickle=True)


def read_calving_timeseries_from_nc(project_dir):

    file_path = os.path.join(project_dir, 'Iceberg Size Distribution.nc')
    ds = nc4.Dataset(file_path)

    glacier_names = list(ds.groups)

    years = [int(year) for year in ds['year'][:]]

    all_timeseries = {}
    for g in range(len(glacier_names)):
        grp = ds.groups[glacier_names[g]]
        all_timeseries[glacier_names[g]]=np.column_stack(
            [
                grp.variables['iceberg_count'][:].T,
                grp.variables['volume'][:]
            ]
        )
    ds.close()

    return all_timeseries, years

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def find_closest_wet_cells(names, locations, fluxes, XC, YC, Depth):

    # fortran starts counting from 1
    rows = np.arange(1,np.shape(XC)[0]+1)
    cols = np.arange(1,np.shape(YC)[1]+1)
    Cols, Rows = np.meshgrid(cols, rows)

    grid_points = np.column_stack([Cols.ravel(), Rows.ravel()])
    model_points = np.column_stack([XC.ravel(), YC.ravel()])
    model_depths = Depth.ravel()

    model_points = model_points[model_depths>100, :]
    grid_points = grid_points[model_depths>100, :]
    model_depths = model_depths[model_depths>100]

    bbox = np.vstack([np.column_stack([XC[0,:], YC[0,:]]),
                      np.column_stack([XC[:,-1], YC[:,-1]]),
                      np.flipud(np.column_stack([XC[-1,:], YC[-1,:]])),
                      np.flipud(np.column_stack([XC[:,0], YC[:,0]]))])
    p = mplPath.Path(bbox)

    # plt.plot(bbox[:,0], bbox[:,1])
    # for ll in range(len(locations)):
    #     plt.plot(locations[ll][0],locations[ll][1],'go')
    # plt.show()

    output_grid = np.zeros((len(locations),6))
    counter = 1
    for ll in range(len(locations)):
        location = locations[ll]
        if p.contains_point((location[0], location[1])):
            distances = great_circle_distance(location[0], location[1], model_points[:,0], model_points[:, 1])
            index = np.argmin(distances)

            # # if overwriting (not allowing multiple glaciers to calve into the same spot)
            # # then use this code
            # if ~np.any(output_grid[:,-1]==index):
            #     row_index = counter-1
            #     output_grid[row_index, 0] = counter
            #     output_grid[row_index, 1] = grid_points[index, 0]
            #     output_grid[row_index, 2] = grid_points[index, 1]
            #     output_grid[row_index, 5] = index
            #     counter += 1

            # if multuple glaciers can calve into the same spot
            row_index = counter - 1
            output_grid[row_index, 0] = counter
            output_grid[row_index, 1] = grid_points[index, 0]
            output_grid[row_index, 2] = grid_points[index, 1]
            output_grid[row_index, 4] = fluxes[ll]/np.sum(fluxes)

            if fluxes[ll]>10:
                print('    '+names[ll]+', '+str(locations[ll])+', '+str(fluxes[ll])+','+str(fluxes[ll]/np.sum(fluxes)))
                print('    Location: '+str(grid_points[index, 0])+','+str(grid_points[index, 1]))

            counter += 1


    output_grid = output_grid[output_grid[:,0]!=0,:]
    output_grid = output_grid[:,:-1]
    # output_grid[:,-1] = 1/np.shape(output_grid)[0]

    return(output_grid)

def look_up_nearest_coastline_point(x, y, coastline_points):
    min_loc = (-1, -1, 100000)
    for row in coastline_points:
        distance = get_distance(row[0], row[1], x, y)
        if distance < min_loc[2]:
            min_loc = (row[0], row[1], distance)
    return int(min_loc[0]), int(min_loc[1])

def create_calving_location_grid(glacier_location_list, coastline_points):

    output_grid = np.zeros((len(glacier_location_list),3))
    counter = 1

    for glacier in range(len(glacier_location_list)):

        row, col = look_up_nearest_coastline_point(glacier_location_list[glacier][0], glacier_location_list[glacier][1], coastline_points)
        # print(glacier+' ('+str(glacier)+')')

        # for i in range(len(rows)):
        # multiple glaciers can calve into the same spot
        row_index = counter - 1
        output_grid[row_index, 0] = counter
        output_grid[row_index, 1] = row+1 # col in fortran counting (!)
        output_grid[row_index, 2] = col+1 # row in fortran counting (!)
        counter += 1

    output_grid = output_grid[output_grid[:,0]!=0,:]

    return(output_grid)

def write_calving_location_file(output_path, output_table):

    output = ''

    for row in range(np.shape(output_table)[0]):
        for col in range(np.shape(output_table)[1]):
            # output+='{:<4.1f}'.format(output_table[row,col])
            if col==0:
                output+=str(output_table[row,col]).rjust(6,' ')
            if col==1:
                output+=str(output_table[row,col]).rjust(7,' ')
            if col==2:
                output+=str(output_table[row,col]).rjust(7,' ')+'   '
        if row<np.shape(output_table)[0]-1:
            output+='\n'

    f = open(os.path.join(output_path, 'calving_locations.txt'),'w')
    f.write(output)
    f.close()

# def write_calving_reference_file(config_dir, glacier_location_dict):
#
#     output = 'Glacier,Row,Col'
#     for glacier in glacier_location_dict:
#         rows = glacier_location_dict[glacier][0]
#         cols = glacier_location_dict[glacier][1]
#         for i in range(len(rows)):
#             output+='\n'+glacier+','+str(rows[i])+','+str(cols[i])
#
#     output_file = os.path.join(config_dir,'L2',model_name,'input','Calving Location Reference File.csv')
#     f = open(output_file,'w')
#     f.write(output)
#     f.close()
#     a=1

def compute_size_distribution(total_volume):

    # set a fixed number of icebergs
    N = 1000

    # compute a normal distribution in log space
    # where the mean is the tuning nob
    # find the mean such that the total volume
    # of the icebergs, based on the size categories, is equal to the flux
    log_v = np.linspace(1, 10, 100)
    v = 10**log_v
    sigma = 2

    mus = np.arange(0,10,0.001)
    err = 1e22
    mu = 1
    for mu_test in mus:
        count = N * (1 / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(-1 * ((log_v - mu_test) ** 2) / (2 * sigma ** 2))
        count = np.round(count).astype(int)
        volume_check = np.sum(v * count)
        err_check = np.abs(volume_check-total_volume)
        if err_check<err:
            err = err_check
            mu = mu_test

    # print(mu)
    count = N*(1/np.sqrt(2*np.pi*sigma**2))*np.exp(-1*((log_v-mu)**2)/(2*sigma**2))
    count = np.round(count).astype(int)

    # print('total number of bergs',np.sum(count))
    volume_check = np.sum(v * count)
    # print('mu',mu)
    # print('        - Volume check: ',total_volume, volume_check,total_volume/volume_check)
    # # print(total_volume/volume_check)
    # max_size_box = np.sum(count!=0)
    # print('max size: '+str(10**log_v[max_size_box])+' m^3')
    # print('max size: ' + str((10 ** log_v[max_size_box])**0.333) + ' m')
    # plt.plot(log_v,count)
    # plt.show()

    return(v, count, mu, sigma)

def date_to_iter_number(date,seconds_per_iter = 60):

    total_seconds = (date-datetime(1992,1,1)).total_seconds()
    iter_number = total_seconds/seconds_per_iter
    # print(iter_number)
    return(iter_number)

def generate_timeseries(years, v, count, N=10000):
    full_times = np.zeros((N*len(years),))
    full_widths = np.zeros((N*len(years),))
    full_thicknesses = np.zeros((N*len(years),))
    full_volumes = np.zeros((N*len(years),))

    start_time = (datetime(years[0],1,1) - datetime(1992,1,1)).total_seconds()
    for year in years:

        if year%4==0:
            end_time = start_time + 366*86400
        else:
            end_time = start_time + 365*86400

        # get a list of all the volumes and the shuffle them
        volumes = np.zeros((N,))
        counted = 0
        for i in range(len(v)):
            volumes[counted:counted+int(count[i])] = v[i]
            counted+=int(count[i])
        np.random.shuffle(volumes)

        # sample the volumes into an estimate of widths and thicknesses
        widths = np.zeros((N,))
        thicknesses = np.zeros((N,))
        removed_for_depth = 0
        for j in range(len(volumes)):
            if volumes[j]!=0:
                cube_width = volumes[j]**0.33

                w = np.random.normal(loc=cube_width,scale=0.2*cube_width)
                t = volumes[j]/(w*1.62*w)
                dep = volumes[j] / (widths[j] * thicknesses[j])
                # rotate if necessary
                if t>w:
                    tmp = np.copy(t)
                    t = np.copy(w)
                    w = tmp

                if t<dep-10:
                    widths[j] = w
                    thicknesses[j]= t
                else:
                    removed_for_depth+=1

        # get calving times
        times = np.round(np.random.uniform(start_time,end_time,(N,)))

        # sort everything by the times
        indices = np.argsort(times)
        times = times[indices]
        widths = widths[indices]
        thicknesses = thicknesses[indices]
        volumes = volumes[indices]

        # remove tiny icebergs
        indices_small = np.logical_and(widths<10, thicknesses<10)
        times[indices_small] = 0
        widths[indices_small] = 0
        thicknesses[indices_small] = 0
        volumes[indices_small] = 0

        # put everything into the big matrix
        full_times[N * (year - years[0]):N * (year - years[0] + 1)] = times
        full_widths[N * (year - years[0]):N * (year - years[0] + 1)] = widths
        full_thicknesses[N * (year - years[0]):N * (year - years[0] + 1)] = thicknesses
        full_volumes[N * (year - years[0]):N * (year - years[0] + 1)] = volumes

        if year%4==0:
            start_time += 366*86400
        else:
            start_time += 365*86400

    return np.column_stack([full_times,full_widths,full_thicknesses,full_volumes])

def create_calving_schedules(output_path, n_bergs_per_year, calving_timeseries, outlet_id):

    np.random.seed(17)

    calving_schedules = np.zeros((int(sum(n_bergs_per_year))*len(years),4))
    for t in range(np.shape(calving_timeseries)[0]):
        dec_yr = calving_timeseries[t,0]
        width = calving_timeseries[t,1]
        length = calving_timeseries[t,2]
        thickness = calving_timeseries[t,3]

        # model_time = (datetime(int(np.floor(dec_yr)),1,1) - datetime(1992,1,1)).total_seconds()
        # if int(np.floor(dec_yr)) %4 == 0:
        #     model_time += (dec_yr-int(np.floor(dec_yr)))*366*24*60*60
        # else:
        #     model_time += (dec_yr - int(np.floor(dec_yr))) * 366 * 24 * 60 * 60
        # date_check = datetime(1992,1,1) + timedelta(seconds=model_time)
        model_time = dec_yr

        calving_schedules[t][0] = model_time
        calving_schedules[t][1] = width
        calving_schedules[t][2] = length
        calving_schedules[t][3] = thickness

    os.makedirs(os.path.join(output_path, 'calving_schedules'), exist_ok=True)
    print('    - Writing schedule for outlet '+str(outlet_id))
    output_schedules = np.array(calving_schedules).T
    output_schedules.ravel('C').astype('>f8').tofile(
        os.path.join(output_path, 'calving_schedules', 'calving_schedule_' + '{:03d}'.format(outlet_id)))

def plot_time_series(data_frame, diff=False):
    fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(20,10))
    if diff:
        for i in range(3):
            col = i + 3
            plot_df = data_frame[data_frame.columns[col]]
            axs[i].plot(plot_df)
            axs[i].set_title(data_frame.columns[col])
    else:
        for i in range(3):
            col = i
            plot_df = data_frame[data_frame.columns[col]]
            axs[i].plot(plot_df)
            axs[i].set_title(data_frame.columns[col])
    plt.show()
    plt.close(fig)


input_folder = '/Volumes/T7/IceBerg/Data/Glacier/'

project_dir = '/Volumes/T7/IceBerg'

print('    - Reading in the model grid')
grid_path = os.path.join('/Users/eddie/IcebergModeling/Input/greenland.mitgrid')

if not os.path.exists(os.path.join(input_folder, 'calving_locations.txt')):
    print('    - Reading in the model ice front locations')
    glacier_location_list = read_glacier_loc_from_loc_file(os.path.join(input_folder, 'gate_loc.npy'))
    converted_glacier_locations = mapping_to_row_and_col(glacier_location_list, grid_path, 360, 180)
    coastline_points = np.load(os.path.join(input_folder, 'coast_line_points.npy'), allow_pickle=True)
    print('    - Writing calving locations file for model')
    output_table = create_calving_location_grid(converted_glacier_locations, coastline_points)
    # print(output_table)
    print('         - Identified '+str(np.shape(output_table)[0])+' calving locations')
    write_calving_location_file(input_folder, output_table)

# print('    - Writing calving reference file (for plots and analysis)')
# output_file = ''
# write_calving_reference_file(config_dir, glacier_location_dict)
#
# # plt.plot(output_table[:,1], output_table[:,2],'k.')
# # plt.show()
#
print('    - Reading in calving timeseries')
all_calving_timeseries, years = read_calving_timeseries_from_nc(input_folder)
if years == []:
    print('    - No years data using default years (2016-2020)')
    years = range(2016, 2021, 1)

glacier_names = list(all_calving_timeseries.keys())
print(f'    - Identified {len(glacier_names)} calving locations')

for i, glacier_name in enumerate(glacier_names):
    print(f'    - Generating calving schedule for {glacier_name}')
    iceberg_count = all_calving_timeseries[glacier_name][0]
    calving_timeseries = generate_timeseries(
        years=years,
        v=all_calving_timeseries[glacier_name][1],
        count=iceberg_count
    )
    print(f'    - Generated calving schedule for {glacier_name}, interpolating to schedule file')
    create_calving_schedules(
        input_folder,
        n_bergs_per_year=iceberg_count,
        calving_timeseries=calving_timeseries,
        # n_glaciers=len(glacier_names),
        outlet_id=i+1
    )













