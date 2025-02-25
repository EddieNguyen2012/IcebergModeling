import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import moviepy.video.io.ImageSequenceClip
import netCDF4 as nc4
import cmocean
from datetime import datetime, timedelta

def extract_time(filename):
    # Extract the date part of the filename
    checkpoint = filename[:-4]  # Extracting the date portion
    # Convert to datetime object
    return int(checkpoint)


def get_full_path_list(file_type):
    frames_folder = os.path.join(file_type)
    files = [f for f in os.listdir(frames_folder) if not f.startswith("._")]
    sorted_files = sorted(files, key=extract_time)
    full_paths = [os.path.join(frames_folder, filename) for filename in sorted_files]
    return full_paths


def create_movie(frame_folder, fps):
    full_paths = get_full_path_list(frame_folder)
    # use the ImageSequenceClip module to set up the clip
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(full_paths, fps=fps)

    # write the video to a file
    clip.write_videofile(os.path.join(frame_folder + '.mp4'), fps=fps)

##################################################################
# Declare constant for functions
n_rows = 360
n_cols = 180
data_path = '/Volumes/T7/output/run/diags/daily_snapshot/daily_snapshot'
years = ['1992']

file_types = [file_type for file_type in os.listdir(data_path) if not file_type.startswith('._')]

variable_metadata = {
    "EtaN": {"vmin": -2.9, "vmax": 2.3, "cmap": "viridis"},
    # "SIarea": {"vmin": 0, "vmax": 1, "cmap": "cmo.ice"},  # Sea ice concentration (0-1)
    # "SIheff": {"vmin": 0, "vmax": 5, "cmap": "cmo.ice"},  # Sea ice thickness (m)
    # "SIhsnow": {"vmin": 0, "vmax": 2, "cmap": "cmo.ice"},  # Snow depth on sea ice (m)
    # "SIvice": {"vmin": -0.5, "vmax": 0.5, "cmap": "cmo.balance"},  # Sea ice velocity
    # "SIuice": {"vmin": -0.5, "vmax": 0.5, "cmap": "cmo.balance"},  # Sea ice u-component velocity
    "Salt_surf": {"vmin": 30, "vmax": 40, "cmap": "cmo.haline"},  # Surface salinity (PSU)
    "Theta_surf": {"vmin": -2, "vmax": 30, "cmap": "cmo.thermal"}, # Sea surface temperature (°C)
    "Salt_AW": {"vmin": 30, "vmax": 40, "cmap": "cmo.haline"},  # Atlantic Water salinity (PSU)
    "Theta_AW": {"vmin": -2, "vmax": 30, "cmap": "cmo.thermal"},  # Atlantic Water temperature (°C)
}


for day in range(61):
    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(2, 3, wspace=0.8, hspace=0.1,
                  left=0.11, right=0.9, top=0.95, bottom=0.05
                  )
    variable_list = list(variable_metadata.keys())
    for i in range(len(variable_metadata)):
        file_type = variable_list[i]
        file = os.listdir(os.path.join(data_path, file_type))[0]
        file_path = os.path.join(data_path, file_type, file)
        ds = nc4.Dataset(file_path)

        grid = ds.variables[file_type][:, :, :]
        plot_grid = np.ma.masked_where(grid == 0.0, grid)

        ax1 = fig.add_subplot(gs[i])
        C = plt.pcolormesh(plot_grid[day, :, :],
                           vmin=variable_metadata[file_type]['vmin'],
                           vmax=variable_metadata[file_type]['vmax'],
                           cmap=variable_metadata[file_type]['cmap'])

        plt.colorbar(C, label=file_type, fraction=0.026)
    plt.savefig(os.path.join('/Users/eddie//IcebergModeling', 'Output', 'Plots', f'{day}.png'))
    plt.close(fig)

create_movie('/Users/eddie/IcebergModeling/Output/Plots', 2)
