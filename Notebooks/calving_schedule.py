import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
    plt.savefig('/Users/eddie/IcebergModeling/Output/Plots/calving_schedule/output.png')
    return fig


# Generate a time series index for one year
original_random = np.random.get_state()


# Define seasonal volatility and glacier thickness bounds (hypothetical values in meters)
def get_glacier_volatility(day_of_year):
    if day_of_year <= 90:  # Q1 (Winter: Ice accumulation, stable)
        return -0.1, 1
    elif day_of_year <= 180:  # Q2 (Spring: Moderate melt begins)
        return 0.3, 10
    elif day_of_year <= 270:  # Q3 (Summer: Heavy ice melt, high fluctuations)
        return -0.1, 10
    else:  # Q4 (Fall: Cooling, ice begins reforming)
        return -0.7, 1


def string_to_ascii(s):
    return sum([ord(char) for char in s])


def convert_time_to_seconds(time):
    reference_date = pd.to_datetime('1992-01-15')
    return pd.to_timedelta(time - reference_date).total_seconds()


# Generate time series with controlled range and volatility
def get_glacier_size(year):
    date_range = pd.date_range(start=str(year) + "-01-15", periods=365, freq="D")
    date_range_in_seconds = convert_time_to_seconds(date_range)
    n = len(date_range)
    dimensions = []
    for dimension in ['Width', 'Length', 'Thickness']:
        np.random.seed(year + string_to_ascii(dimension))
        ts_values = [np.random.randint(200, 225, 1)[0]]  # Start with an initial glacier thickness of 60 meters
        np.random.set_state(original_random)
        for i in range(1, n):
            day_of_year = date_range[i].dayofyear
            roc, volatility = get_glacier_volatility(day_of_year)

            noise = np.random.normal(0, volatility)  # Scale noise
            new_value = roc + ts_values[i - 1] + noise  # Accumulate changes

            # Enforce the range
            # new_value = np.clip(new_value, ts_values[0] - )
            new_value = np.clip(new_value, 20, 500)
            ts_values.append(new_value)

        # Create Pandas Series
        dimensions.append(pd.Series(ts_values, index=date_range_in_seconds, name=dimension))

    dimensions = pd.concat(dimensions, axis=1)
    return dimensions

#------------------main program----------------------


output_array = np.array(np.zeros((100000, 4)))
multiple_timeseries = []
for i in range(1992, 1997):
    multiple_timeseries.append(get_glacier_size(i))
final_df = pd.concat(multiple_timeseries)
final_array = final_df.reset_index().to_numpy()
index = 3
print(output_array)
for i in range(len(final_df)):
    output_array[index, 0] = final_array[i, 0]
    output_array[index, 1] = final_array[i, 1]
    output_array[index, 2] = final_array[i, 2]
    output_array[index, 3] = final_array[i, 3]

    index += np.random.randint(1, 10, 1)

output_array = output_array.T
print(output_array)
for i in range(1, 6):
    output_array.ravel('C').astype('>f8').tofile(f'calving_schedule_00{i}')

plot_time_series(final_df)