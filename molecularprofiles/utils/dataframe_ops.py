import pandas as pd
import numpy as np
import colorlover as cl
import plotly.plotly as py
import plotly.io as pio
import plotly.graph_objs as go
from molecularprofiles.utils.observatory import *


def select_dataframe_epoch(df, epoch_text):
    epoch = get_epoch(epoch_text)
    new_df = df[df.month.isin(epoch)]
    return new_df

    # if epoch_text == 'winter':
    #     condition = (df.month == epoch[0]) | (df.month == epoch[1]) | (df.month == epoch[2]) | (df.month == epoch[3])
    #     new_df = df[condition]
    #
    # elif epoch_text == 'summer':
    #     condition = (df.month == epoch[0]) | (df.month == epoch[1]) | (df.month == epoch[2])
    #     new_df = df[condition]
    #
    # elif epoch_text == 'intermediate':
    #     condition = (df.month == epoch[0]) | (df.month == epoch[1]) | (df.month == epoch[2]) | (df.month == epoch[3]) | \
    #                 (df.month == epoch[4])
    #     new_df = df[condition]
    # return new_df


def select_dataframe_by_year(df, years):
    new_df = df[df.year.isin(years)]
    return new_df


def select_dataframe_by_month(df, months):
    new_df = df[df.month.isin(months)]
    return new_df


def select_dataframe_by_hour(df, hours):
    new_df = df[df.hour.isin(hours)]
    return new_df


def create_wind_speed_dataframe(df, normalized=False):
    """
    Function to create a wind speed dataframe in order to plot it afterwards as a wind rose plot

    :param df: dataframe containing wind direction information
    :param normalized (optional):
    :return: df_winds
    """
    wd_centre_bins = np.arange(7.5,360, 15)
    ws_hist = []
    for d in wd_centre_bins:
        ws_hist.append(np.histogram(df.wind_speed[(df.wind_direction >= d-7.5) & (df.wind_direction < d + 7.5)],
                                    bins=[0,5,10,20,30,40,50,100])[0])

    df_winds = pd.DataFrame(columns=['wind_direction', '0-5', '5-10', '10-20', '20-30', '30-40', '40-50', '> 50'])
    ws_new_list = []
    for j in range(len(ws_hist[0])):
        li = []
        for i in range(len(ws_hist)):
            li.append(ws_hist[i][j])
        ws_new_list.append(li)

    for i,j in zip(df_winds.keys()[1:], range(len(ws_new_list))):
        df_winds[i] = ws_new_list[j]
    if normalized:
        df_winds_normalized = df_winds.div(df_winds.sum(axis=1), axis=0)
        df_winds_normalized['wind_direction'] = wd_centre_bins
        return df_winds_normalized
    else:
        df_winds['wind_direction'] = wd_centre_bins
        return df_winds


def plot_wind_rose(df, name_tag='my_wind_rose'):
    """
    Function to plot a wind rose graphic from a pandas dataframe object that contains wind direction
    information
    :param df: a pandas dataframe created with the function create_wind_speed_dataframe
    :param name_tag (optional):
    :return:
    """
    data = []
    counter = 0
    for col in df.columns:
        if col != 'wind_direction':
            data.append(go.Area(t=df['wind_direction'], r=df[col],
                                marker=dict(color=cl.scales['9']['seq']['YlGnBu'][counter]), name=col+' m/s'))
            counter+=1
    print(data)

    fig = go.Figure(data=data, layout=go.Layout(orientation=270., barmode='stack'))
    pio.write_image(fig, name_tag+'_wind_speed_rose.pdf')
    py.offline.plot(fig)


def compute_averages_std_simple(input_array):
    """
    This function computes the average, standard deviation and peak to peak (plus and minus)
    for an input 1D array

    Input:
        1-D input_array (array-like)

    Output:
        average, stand_dev, peak_to_peak_p, peak_to_peak_m (array-like)
    """
    import numpy as np
    from statsmodels import robust

    average = np.average(input_array)
    stand_dev = robust.mad(input_array)
    peak_to_peak_p = np.max(input_array) - np.average(input_array)
    peak_to_peak_m = np.average(input_array) - np.min(input_array)
    return average, stand_dev, peak_to_peak_p, peak_to_peak_m


def avg_std_dataframe(group, param):
    """

    :param group: dataframe grouped by a certain parameter
    :param param: the parameter by which the dataframe is grouped
    :return:
        avg: the mean value for each grouped level
        std: the standard deviation for each grouped level
        mad: the mean absolute deviation for each group level
        p2p_p: the peak-to-peak maximum value for each grouped level
        p2p_m: the peak-to-peak minimum value for each grouped level
    """

    avg = group[param].mean()
    std = group[param].std()
    mad = group[param].mad()
    med = group[param].median()
    p2p_p = np.max(group[param]) - avg
    p2p_m = avg - np.min(group[param])

    return avg, std, mad, p2p_p, p2p_m


def compute_averages_std(input_array):
    """
    This function computes the average, standard deviation and peak to peak (plus and minus)
    for a multidimensional input array

    Input: input_array (array-like)

    Output: average, stand_dev, peak_to_peak_p, peak_to_peak_m (all array-like)
    """
    import numpy as np
    from statsmodels import robust
    average = []
    stand_dev = []
    peak_to_peak_p = []
    peak_to_peak_m = []

    for i in np.arange(len(input_array[0])):
        average.append(np.average(input_array[:, i]))
        stand_dev.append(robust.mad(input_array[:, i]))
        peak_to_peak_p.append(np.max(input_array[:, i] - np.average(input_array[:, i])))
        peak_to_peak_m.append(np.average(input_array[:, i]) - np.min(input_array[:, i]))

    average = np.asarray(average)
    stand_dev = np.asarray(stand_dev)
    peak_to_peak_p = np.asarray(peak_to_peak_p)
    peak_to_peak_m = np.asarray(peak_to_peak_m)
    return average, stand_dev, peak_to_peak_p, peak_to_peak_m
