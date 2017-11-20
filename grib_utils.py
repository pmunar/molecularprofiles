from builtins import str

def angle2radians(deg, min, sec):
    """
    This function computes the value, in radians, of an angle expressed in degrees, minutes, seconds
    Input: (integers) degrees, minutes, seconds
    Output: (float) angle in radians
    """
    import math
    if deg > 0.:
        angle = deg + (min / 60.) + (sec / 3600.)
    else:
        angle = deg - (min / 60.) - (sec / 3600.)
    return math.radians(angle)

def get_observatory_coordinates(observatory):
    if observatory == 'north':
        latitude = angle2radians(28, 45, 42.462)
        longitude = angle2radians(-17, 53, 26.525)
    elif observatory == 'south':
        latitude = angle2radians(-24, 40, 24.8448)
        longitude = angle2radians(-70, 18, 58.4712)
    else:
        print('wrong observatory!')
        exit()

    return latitude, longitude

def GetAltitudeFromGeopotentialHeight(geop_height, observatory):
    import numpy as np
    latitude = get_observatory_coordinates(observatory)[0]
    geop_heightkm = geop_height / 1000. / 9.80665  # dividing by the acceleration of gravity on Earth
    cos2lat = np.cos(2 * latitude)
    # convert from geopotential height to geometric altitude:
    z = (1. + 0.002644 * cos2lat) * geop_heightkm + (1 + 0.0089 * cos2lat) * geop_heightkm * geop_heightkm / 6245.
    # convert Z to meter
    return 1.0E3 * z  # This is fGeoidOffset


def date2mjd(yyyy, mm, dd, hh, min=0, ss=0):
    """
    This function computes the mjd value corresponding to an input date, in the year, month, day, hour format

    Input:
        (integers) year, month, day and hour.
        Optionals: minutes and seconds. If left blank, they are assumed equal 0
    
    Output:
        (float) mjd
    """
    from astropy.time import Time
    y = str(int(yyyy))
    m = '{0:02}'.format(int(mm))
    d = '{0:02}'.format(int(dd))
    h = '{0:02}'.format(int(hh))
    mi = '{0:02}'.format(int(min))
    s = '{0:02}'.format(int(ss))
    t = Time(y + '-' + m + '-' + d + 'T' + h + ':' + mi + ':' + s, format='isot', scale='utc')
    return round(t.mjd, 2)


def mjd2date(mjd):
    """
    This function computes the date from an input mjd day

     Input: mjd (float)

    Output: year, month, day, hour (integers) corresponding to the mjd input day
    """
    from astropy.time import Time
    time = Time(mjd, format='mjd')
    return int(time.datetime.year), int(time.datetime.month), int(time.datetime.day), int(time.datetime.hour)

def compute_averages_std_simple(input_array):
    """
    This function computes the average, standard deviation and peak to peak (plus and minus) for an input array

    Input:
        input_array (array-like)

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


def compute_averages_std(input_array):
    """
    This function computes the average, standard deviation and peak to peak (plus and minus)
    for a multidimensional input array

    Input: input_array (array-like)

    Output: average, stand_dev, peak_to_peak_p, peak_to_peak_m (array-like)
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

def get_winter_months():
    import numpy as np
    return np.array([1, 2, 3, 4, 12])


def get_summer_months():
    import numpy as np
    return np.array([7, 8, 9])


def get_all_months():
    import numpy as np
    return np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])


def get_intermediate_months():
    import numpy as np
    return np.array([5, 6, 10, 11])


def get_epoch(epoch):
    valid_epochs = {'winter' : get_winter_months(),
                  'summer' : get_summer_months(),
                  'all' : get_all_months(),
                  'intermediate' : get_intermediate_months()}
    return valid_epochs[epoch]

# def print_grib_to_txt(file_name, ):
#     table_file = open(file_name[:-5] + '.txt', 'w')
#     print('# Date, year, month, day, hour, Plevel, T_average, geopotential_average, height, n', file=table_file)
#     print('creating the txt file containing the selected data...')
#     for j in np.arange(len(temperature)):
#         print(temperature[j].dataDate, temperature[j].year, temperature[j].month, temperature[j].day,
#               temperature[j].hour, temperature[j].level, temperature_gp, geopotential_gp, h, density, file=table_file)
#     table_file.close()

def get_closest_gridpoint(lat, lon):
    """

    :param lat: float
    :param lon: float
    :return: nearest ecmwf grid point
    """
    step = 0.75  # grid step in degrees
    lons_grid = np.zeros(int(360 / step) + 1)
    lats_grid = np.zeros(int(180 / step) + 1)

    # getting the grid points for ecmwf data:
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


def readgribfile2text(file_name, year, observatory):
    """
    This function opens a grib file, selects the Temperature and Geopotential parameters,
    and finally creates a txt file where these parameters, together with date, year, month,
    day, hour, pressure level, real height and density, are written.

    Input: file_name (string)

    Output: a txt file with the exact name as the input file name, but with .txt as extension
    """
    import pygrib as pg
    import numpy as np

    Ns = 2.546899E19  # [cm^-3] molecular number density for standard air conditions
    ps = 1013.25      # [hPa]   standard air pressure
    Ts = 288.15       # [K]     standard air temperature
    Hs = 9500.        # [km]    density scale hight for La Palma Winter

    latitude_obs, longitude_obs = get_observatory_coordinates(observatory)

    lat_gridpoint, lon_gridpoint = get_closest_gridpoint(latitude_obs, longitude_obs)

    print('opening grib file (this might take a while...)')
    grb = pg.open(file_name)
    print('selecting temperature parameter...')
    temperature = grb.select(name='Temperature')

    print('selecting geopotential parameter...')
    geop_height = grb.select(name='Geopotential')

    # We create the table file and fill it with the information stored in the above variables, plus the height
    # and density computed form them.
    table_file = open(file_name[:-5] + '.txt', 'w')
    print('# Date, year, month, day, hour, Plevel, T_average, geopotential_average, height, n', file=table_file)

    print('creating the txt file containing the selected data...')

    for j in np.arange(len(temperature)):
        if type(temperature[j].values) == float:
            geopotential_gp = geop_height[j].values #[geop_height[j].year == year]
            temperature_gp = temperature[j].values #[temperature[j].year == year]
        else:
            geopotential_gp = np.float(geop_height[j].values[(geop_height[j].data()[1] == lat_gridpoint) & (
            geop_height[j].data()[2] == lon_gridpoint) & (geop_height[j].year == year)])
            temperature_gp = np.float(temperature[j].values[(temperature[j].data()[1] == lat_gridpoint) & (
            temperature[j].data()[2] == lon_gridpoint) & (temperature[j].year == year)])
        h = GetAltitudeFromGeopotentialHeight(geopotential_gp, observatory)
        density = Ns * temperature[j].level / ps * Ts / temperature_gp
        print(temperature[j].dataDate, temperature[j].year, temperature[j].month, temperature[j].day,
              temperature[j].hour, temperature[j].level, temperature_gp, geopotential_gp, h, density, file=table_file)
    table_file.close()
