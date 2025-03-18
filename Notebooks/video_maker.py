import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import moviepy.video.io.ImageSequenceClip
import netCDF4 as nc4
import cmocean
from datetime import datetime, timedelta


def extract_time(filename):
    # Extract the date part of the filename
    checkpoint = filename[-12:-4]  # Extracting the date portion
    # Convert to datetime object
    return int(checkpoint)


def get_full_path_list(input_dir):
    files = [f for f in os.listdir(input_dir) if not f.startswith("._")]
    sorted_files = sorted(files, key=extract_time)
    full_paths = [os.path.join(input_dir, filename) for filename in sorted_files]
    return full_paths


def create_movie_clip(frame_folder, fps):
    full_paths = get_full_path_list(frame_folder)
    # use the ImageSequenceClip module to set up the clip
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(full_paths, fps=fps)
    return clip
    # write the video to a file


def date_range_extraction(sample_list):
    time_list = set()
    for sample in sample_list:
        year = int(sample.split('_')[-1][:4])
        month = int(sample.split('_')[-1][4:6])
        tmp_timeframe = TimeFrame(year=year, month=month)
        time_list.add(tmp_timeframe)
    return list(sorted(time_list))


class TimeFrame:
    def __init__(self, year: int, month: int):
        self.year = year
        self.month = month

    def __lt__(self, other):
        return (self.year, self.month) < (other.year, other.month)


def date_range_extraction_with_range(sample_list, start_date: datetime, end_date: datetime):
    time = set()
    for sample in sample_list:
        year = int(sample.split('_')[-1][:4])
        month = int(sample.split('_')[-1][4:6])
        tmp_date = TimeFrame(year=year, month=month)
        if TimeFrame(year=start_date.year, month=start_date.month) <= tmp_date <= TimeFrame(year=end_date.year,
                                                                                            month=end_date.month):
            time.add(tmp_date)
    return list(sorted(time))


class Variable:
    def __init__(self, variable: str):
        self.variable = variable
        self.file_count = 0

    def count_file(self, input_dir):
        try:
            self.file_count = len(os.listdir(os.path.join(input_dir, self.variable)))
        except FileNotFoundError:
            self.file_count = 0

    def __lt__(self, other):
        return self.file_count < other.file_count

    def __str__(self):
        return f'var: {self.variable}, count: {self.file_count}'


class MovieConfiguration:
    variable_metadata = {
        "EtaN": {"vmin": -2.9, "vmax": 2.3, "cmap": "viridis", "unit": "$^{\circ}$C"},
        "SIarea": {"vmin": 0, "vmax": 1, "cmap": "cmo.ice", "unit": ""},  # Sea ice concentration (0-1)
        "SIheff": {"vmin": 0, "vmax": 5, "cmap": "cmo.ice", "unit": "m"},  # Sea ice thickness (m)
        "SIhsnow": {"vmin": 0, "vmax": 2, "cmap": "cmo.ice", "unit": "m"},  # Snow depth on sea ice (m)
        "SIvice": {"vmin": -0.5, "vmax": 0.5, "cmap": "cmo.balance", "unit": "m/s"},  # Sea ice velocity
        "SIuice": {"vmin": -0.5, "vmax": 0.5, "cmap": "cmo.balance", "unit": "m/s"},  # Sea ice u-component velocity
        "Salt_surf": {"vmin": 30, "vmax": 40, "cmap": "cmo.haline", "unit": "PSU"},  # Surface salinity (PSU)
        "Theta_surf": {"vmin": -2, "vmax": 30, "cmap": "cmo.thermal", "unit": "$^{\circ}$C"},
        # Sea surface temperature (°C)
        "Salt_AW": {"vmin": 30, "vmax": 40, "cmap": "cmo.haline", "unit": "PSU"},  # Atlantic Water salinity (PSU)
        "Theta_AW": {"vmin": -2, "vmax": 30, "cmap": "cmo.thermal", "unit": "$^{\circ}$C"},
        # Atlantic Water temperature (°C)
        "Salt": {"vmin": 30, "vmax": 40, "cmap": "cmo.haline", "unit": "PSU"},  # Surface salinity (PSU)
        "Theta": {"vmin": -2, "vmax": 30, "cmap": "cmo.thermal", "unit": "$^{\circ}$C"},
        "Uvel": {"vmin": -0.5, "vmax": 0.5, "cmap": "cmo.balance", "unit": "m/s"},  # Sea ice velocity
        "Vvel": {"vmin": -0.5, "vmax": 0.5, "cmap": "cmo.balance", "unit": "m/s"},  # Sea ice u-component velocity
    }

    def __init__(self, input_folder: str, output_folder: str, row=360, col=180, plot_width=15, plot_height=12):
        """
        Initialize the configurations of the movie

        :param input_folder: Input folder path (the path to 'diags' folder) :param output_folder: Output folder path
        (the path to the folder that store plots and video) (will create one if it doesn't exist) :param row: The
        number of rows of the bathymetry grid (default: 360) :param col: The number of columns of the bathymetry grid
        (default: 180) :param plot_width: The width of each plot/output movie (default: 15) :param plot_height: The
        height of each plot/output movie (default: 12)
        """
        self.row = row
        self.column = col
        self.plot_width = plot_width
        self.plot_height = plot_height
        self.input_folder = input_folder

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            print(f"Created a new output folder at {output_folder}")

        self.output_folder = output_folder

    def __str__(self):
        return (f'Movie properties: '
                f'row={self.row}, '
                f'column={self.column}, '
                f'plot_width={self.plot_width}, '
                f'plot_height={self.plot_height}, '
                f'input_folder={self.input_folder}, '
                f'output_folder={self}')

    def change_dataset_dimensions(self, new_row, new_col):
        if new_row and new_col:
            self.row = new_row
            self.column = new_col
        else:
            raise ValueError("Missing arguments, please provide a new row count AND a new column count")

    def change_plot_size(self, new_width, new_height):
        if new_width and new_height:
            self.plot_width = new_width
            self.plot_height = new_height
        else:
            raise ValueError("Missing arguments, please provide a new row count AND a new column count")

    def create_plots(self, option: str, start_date: datetime = None, end_date: datetime = None):
        if not option:
            raise ValueError(
                "No plot choice provided, please provide at least one of daily_mean, daily_snapshot, monthly_mean or "
                "monthly_snapshot")

        variables_tmp = set([Variable(folder) for folder in os.listdir(os.path.join(self.input_folder, option)) if
                             (not folder.startswith("._")) and (
                                 os.listdir(os.path.join(self.input_folder, option, folder)))])
        if not variables_tmp:
            raise ValueError("There is no variables in the folder")

        for variable in variables_tmp:
            variable.count_file(os.path.join(self.input_folder, option))

        variables_tmp = sorted(variables_tmp, reverse=True)

        variables = []
        for variable in variables_tmp:
            variables.append(variable.variable)

        if start_date and end_date:
            time_frames = date_range_extraction_with_range(
                [file for file in sorted(os.listdir(os.path.join(self.input_folder, option, variables[0]))) if
                 not file.startswith('._')], start_date, end_date)
        else:
            time_frames = date_range_extraction(
                [file for file in sorted(os.listdir(os.path.join(self.input_folder, option, variables[0]))) if
                 not file.startswith('._')])

        if not os.path.exists(os.path.join(self.output_folder, 'plots', option)):
            os.makedirs(os.path.join(self.output_folder, 'plots', option))
            print(f"Created folder for {option} plots at {os.path.join(self.output_folder, 'plots', option)}")
        for time_frame in time_frames:
            plot_title = option.replace('_', ' ').capitalize()

            data_list = {}
            for variable in variables:
                try:
                    path = os.path.join(self.input_folder, option, variable,
                                        f'{variable}_{time_frame.year:04}{time_frame.month:02}.nc')
                    data_list[variable] = path
                except FileNotFoundError:
                    continue

            # Assume all variables are available in the same date ranges
            tmp_ds = nc4.Dataset(data_list[variables[0]])
            date_count = len(tmp_ds.dimensions['iterations'])
            tmp_ds.close()
            print(f'Creating plot for {time_frame.year:04}/{time_frame.month:02}')
            for date in range(date_count):
                fig = plt.figure(figsize=(15, 12))
                if len(variables) > 1:
                    gs = GridSpec(int(np.ceil(len(variables) / 3)), 3, wspace=0.8, hspace=0.1,
                                  left=0.11, right=0.9, top=0.95, bottom=0.05
                                  )
                else:
                    gs = GridSpec(1, 1, wspace=0.8, hspace=0.1,
                                  left=0.11, right=0.9, top=0.95, bottom=0.05
                                  )
                day = date + 1

                fig.suptitle(f'{plot_title} at {time_frame.month:02}/{day:02}/{time_frame.year:04}')

                # plot
                for i in range(len(variables)):
                    var = variables[i]
                    try:
                        ds = nc4.Dataset(data_list[var])
                    except FileNotFoundError:
                        continue
                    grid = ds.variables[var][:, :, :]
                    plot_grid = np.ma.masked_where(grid == 0.0, grid)

                    ax1 = fig.add_subplot(gs[i])
                    C = plt.pcolormesh(
                        plot_grid[date, :, :],
                        vmin=self.variable_metadata[var]['vmin'],
                        vmax=self.variable_metadata[var]['vmax'],
                        cmap=self.variable_metadata[var]['cmap'])

                    plt.colorbar(C, orientation='vertical',
                                 label=f"{var} \n ({self.variable_metadata[var]['unit']})")
                    ds.close()

                plt.savefig(os.path.join(self.output_folder, 'plots', option,
                                         f'{option}_{time_frame.year:04}{time_frame.month:02}{day:02}.png'))
                plt.close(fig)

        print(f"Successfully created plots for {option}")

    def create_movie(self, option, fps=30):
        if not option:
            raise ValueError(
                "No plot choice provided, please provide at least one of daily_mean, daily_snapshot, monthly_mean or "
                "monthly_snapshot")
        clip = create_movie_clip(os.path.join(self.output_folder, 'plots', option), fps)
        if not os.path.exists(os.path.join(self.output_folder, 'video', option)):
            os.makedirs(os.path.join(self.output_folder, 'video', option))
            print(f"Created folder for {option} videos at {os.path.join(self.output_folder, 'video', option)}")
        clip.write_videofile(os.path.join(self.output_folder, 'video', option, f'{option}_video' + '.mp4'), fps=fps)

################################################################## Main Program functions

#
# def check_input_argument() -> MovieConfiguration:
#     file_choices = []
#     input_dir = ''
#     output_dir = 'movie_output/'
#     parser = argparse.ArgumentParser(prog='video_maker.py', description='Video Maker script for MITgcm model output')
#     parser.add_argument('-dm', '--daily_mean', help='Select Daily Mean outputs', action='store_true')
#     parser.add_argument('-ds', '--daily_snapshot', help='Select Daily Snapshot outputs', action='store_true')
#     parser.add_argument('-mm', '--monthly_mean', help='Select Monthly Mean outputs', action='store_true')
#     parser.add_argument('-ms', '--monthly_snapshot', help='Select Monthly Snapshot outputs', action='store_true')
#     parser.add_argument('-id', '--input_directory', help='Enter input directory (path to "diags" folder)')
#     parser.add_argument('-od', '--output_directory', help='Enter output directory (path to your desired output folder '
#                                                           'to store the movie)')
#     args = parser.parse_args()
#
#     if args.daily_mean:
#         file_choices.append('daily_mean')
#     if args.daily_snapshot:
#         file_choices.append('daily_snapshot')
#     if args.monthly_mean:
#         file_choices.append('monthly_mean')
#     if args.monthly_snapshot:
#         file_choices.append('monthly_snapshot')
#
#     if not any([args.daily_mean, args.daily_snapshot, args.monthly_mean, args.monthly_snapshot]):
#         print("Please provide at least one of --daily_mean, --daily_snapshot, --monthly_mean or --monthly_snapshot")
#
#     if args.input_directory:
#         print("Directory is", args.input_directory)
#     else:
#         print("Please provide input directory with option -d")
#
#     if args.output_directory:
#         output_directory = args.output_directory
#
#     return configuration(input_dir, output_dir)
#
#
# if __name__ == "__main__":
#     check_input_argument()
