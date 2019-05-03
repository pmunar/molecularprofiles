import numpy as np

def ddmmss2deg(deg, min, sec):
    """
    This function computes the value, in degrees, of an angle expressed in degrees, minutes, seconds
    Input: (integers) degrees, minutes, seconds
    Output: (float) angle in degrees
    """

    if deg > 0.:
        angle = deg + (min / 60.) + (sec / 3600.)
    else:
        angle = deg - (min / 60.) - (sec / 3600.)
    return angle



def get_observatory_coordinates(observatory):
    if observatory == 'north':
        latitude = ddmmss2deg(28, 45, 42.462)
        longitude = ddmmss2deg(-17, 53, 26.525)
        return latitude, longitude
    elif observatory == 'south':
        latitude = ddmmss2deg(-24, 40, 24.8448)
        longitude = ddmmss2deg(-70, 18, 58.4712)
        return latitude, longitude
    else:
        print('Unknown observatory!')
        latitude = np.float(input('Observatory latitude (in degrees):'))
        longitude = np.float(input('Observatory longitude (in degrees):'))
        return latitude, longitude


def get_winter_months():
    return [1, 2, 3, 4, 11, 12]


def get_summer_months():
    return [6, 7, 8, 9, 10]


def get_all_months():
    return [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]


def get_intermediate_months():
    return [5, 6, 10, 11]


def get_epoch(epoch):
    valid_epochs = {'winter': get_winter_months(),
                    'summer': get_summer_months(),
                    'intermediate': get_intermediate_months(),
                    'all': get_all_months()}
    return valid_epochs[epoch]


def get_south_winter_months():
    return [5,6,7,8,9,10]


def get_south_summer_months():
    return [1,2,3,4,5,10,11,12]


def get_south_epoch(epoch):
    valid_epochs = {'winter': get_south_winter_months(),
                    'summer': get_south_summer_months(),
                    'all': get_all_months()}
    return valid_epochs[epoch]


def select_new_epochs_dataframe_north(df,epoch_text):

    epoch = get_epoch(epoch_text)

    if epoch_text == 'winter':
        condition = (df.month == epoch[0]) | (df.month == epoch[1]) | (df.month == epoch[2]) | \
                    (df.month == epoch[3]) | ((df.month == epoch[4]) & (df.day > 15)) | (df.month == epoch[5])
        new_df = df[condition]

    elif epoch_text == 'summer':
        condition = ((df.month == epoch[0]) & (df.day > 20)) | (df.month == epoch[1]) | (df.month == epoch[2]) | \
                    (df.month == epoch[3])
        new_df = df[condition]

    elif epoch_text == 'intermediate':
        condition = (df.month == epoch[0]) | ((df.month == epoch[1]) & (df.day <= 20)) | \
                    ((df.month == epoch[2]) & (df.day > 5)) | ((df.month == epoch[3]) & (df.day <= 15))
        new_df = df[condition]

    return new_df


def select_new_epochs_dataframe_south(df,epoch_text):

    epoch = get_south_epoch(epoch_text)

    if epoch_text == 'summer':
        condition = (df.month == epoch[0]) | (df.month == epoch[1]) | (df.month == epoch[2]) | (df.month == epoch[3]) |\
                    ((df.month == epoch[4]) & (df.day < 15)) | ((df.month == epoch[5]) & (df.day > 15)) | \
                    (df.month == epoch[6]) | (df.month == epoch[7])
    elif epoch_text == 'winter':
        condition = ((df.month == epoch[0]) & (df.day > 15)) | (df.month == epoch[1]) | (df.month == epoch[2]) | \
                    (df.month == epoch[3]) | (df.month == epoch[4]) | ((df.month == epoch[5]) & (df.day < 15))
    new_df = df[condition]
    return new_df

def select_new_epochs_dataframe_density_north(df, epoch_text):
    if epoch_text == 'summer':
        condition = df[(df.n_exp > 0.88) & (df.P == 125)]
    elif epoch_text == 'winter':
        condition = df[(df.n_exp < 0.8375) & (df.P == 125)]
    elif epoch_text == 'intermediate':
        condition = df[(df.n_exp < 0.88) & (df.n_exp > 0.8375) & (df.P == 125)]
    new_df = df[df.MJD.isin(condition.MJD)]
    return new_df

def select_new_epochs_dataframe_density_south(df, epoch_text):
    if epoch_text == 'summer':
        condition = df[(df.n_exp > 0.88) & (df.P == 125)]
    elif epoch_text == 'winter':
        condition = df[(df.n_exp < 0.88) & (df.P == 125)]
    new_df = df[df.MJD.isin(condition.MJD)]
    return new_df

def select_epoch(file, epoch_text):
    global months, month
    print("loading and selecting data")
    file = open(file)
    date, year, month, day, hour, mjd, p, T, h, n, U, V, RH = np.loadtxt(file, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                                        10, 11, 12), unpack=True)

    epoch = get_epoch(epoch_text)

    if epoch_text == 'winter':
        condition = (month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) | (month == epoch[3])
        date = date[condition]
        mjd = mjd[condition]
        h = h[condition]
        n = n[condition]
        p = p[condition]
        T = T[condition]
        V = V[condition]
        U = U[condition]
        RH = RH[condition]

    elif epoch_text == 'summer':
        condition = (month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])
        date = date[condition]
        mjd = mjd[condition]
        h = h[condition]
        n = n[condition]
        p = p[condition]
        T = T[condition]
        V = V[condition]
        U = U[condition]
        RH = RH[condition]

    elif epoch_text == 'intermediate':
        condition = (month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) | (month == epoch[3]) | \
                    (month == epoch[4])
        date = date[condition]
        mjd = mjd[condition]
        h = h[condition]
        n = n[condition]
        p = p[condition]
        T = T[condition]
        V = V[condition]
        U = U[condition]
        RH = RH[condition]

    elif epoch_text == 'all':
        date = date
        mjd = mjd
        h = h
        n = n
        p = p
        T = T
        U = U
        V = V
        RH = RH

    return date, year, month, day, hour, mjd, p, h, n, T, U, V, RH


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


def get_closest_gridpoint(lat, lon, gridstep):
    """
    :param lat: float
    :param lon: float
    :param gridstep: float. Step in which the grid is divided (0.75 degrees for ECMWF data
    and 1.0 degrees for GDAS data)
    :return: nearest grid point
    """
    step = gridstep  # grid step in degrees
    lons_grid = np.zeros(int(360 / step) + 1)
    lats_grid = np.zeros(int(180 / step) + 1)

    # getting the grid points:
    for i in range(len(lons_grid)):
        lons_grid[i] = step * i
    for i in range(len(lats_grid)):
        lats_grid[i] = -90 + step * i

    print('Latidude of interest is %5.2f deg' % lat)
    print('Longitude of interest is %5.2f deg' % (lon % 360))
    nearest_lat = find_nearest(lats_grid, lat)
    print('nearest latitude is: %4.1f deg' % nearest_lat )
    nearest_lon = find_nearest(lons_grid, lon % 360)
    print('nearest longitude is: %5.1f deg' % nearest_lon )
    return nearest_lat, nearest_lon


def find_nearest(a, value):  # Function to find the nearest grid position to a given latitude or longitude
    return a[np.abs(a-value).argmin()]
