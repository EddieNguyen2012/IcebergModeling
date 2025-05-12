import argparse

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_time_series(data_frame, schedule_name, diff=False):
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
    plt.suptitle(f'Time series for {schedule_name}')
    plt.savefig(f'{schedule_name}_timeseries_plot.png')
    plt.close(fig)

def convert_seconds_to_time(seconds):
    reference_date = pd.to_datetime('1992-01-15')
    return pd.to_datetime(reference_date + pd.to_timedelta(seconds, unit='s'))

def get_time_series(file_path):
    schedule = np.fromfile(file_path, '>f8')
    schedule = schedule.reshape((4, int(len(schedule) / 4),)).T
    schedule = schedule[np.where(schedule[:, 0] != 0)]
    schedule = pd.DataFrame(schedule, columns=['year', 'width', 'height', 'thickness'])
    schedule['year'] = convert_seconds_to_time(schedule['year'])
    schedule.set_index('year', inplace=True)

    return schedule


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--path", action="store",
                        help="The directory to calving_schedule file.", dest="file_path",
                        type=str, required=True)

    args = parser.parse_args()
    file_path = args.file_path

    schedule_df = get_time_series(file_path)
    plot_time_series(data_frame=schedule_df, schedule_name=str.split(file_path, '/')[-1])