#!/home/pmunar/anaconda3/bin/python3.6

from builtins import str
import numpy as np
import pygrib as pg
from molecularprofiles.utils.meteorological_constants import *
import sys
import gc
import multiprocessing
from multiprocessing import Process
import os

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

def GetAltitudeFromGeopotential(geop_height, observatory):
    """
    Function to compute the real altitude from the geopotential value at a certain coordinates on Earth
    :param geop_height:
    :param observatory: possible values are 'north' or 'south'
    :return: real altitude as fGeoidOffset (in m)
    """
    latitude = np.radians(get_observatory_coordinates(observatory)[0])
    geop_heightkm = geop_height / 1000. / 9.80665  # dividing by the acceleration of gravity on Earth
    cos2lat = np.cos(2 * latitude)
    # convert from geopotential height to geometric altitude:
    z = (1. + 0.002644 * cos2lat) * geop_heightkm + (1 + 0.0089 * cos2lat) * geop_heightkm * geop_heightkm / 6245.
    # convert Z to meter
    return 1.0E3 * z  # This is fGeoidOffset


def GetAltitudeFromGeopotentialHeight(geop, observatory):
    """
    Function to compute the real altitude from the geopotential value at a certain coordinates on Earth
    :param geop_height:
    :param observatory: possible values are 'north' or 'south'
    :return: real altitude as fGeoidOffset (in m)
    """
    latitude = np.radians(get_observatory_coordinates(observatory)[0])
    geop_km = geop / 1000.
    cos2lat = np.cos(2 * latitude)
    # convert from geopotential height to geometric altitude:
    z = (1. + 0.002644 * cos2lat) * geop_km + (1 + 0.0089 * cos2lat) * geop_km * geop_km / 6245.
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
        p2p_p: the peak-to-peak maximum value for each grouped level
        p2p_m: the peak-to-peak minimum value for each grouped level
    """

    avg = group[param].mean()
    std = group[param].std()
    p2p_p = np.max(group[param]) - avg
    p2p_m = avg - np.min(group[param])

    return avg, std, p2p_p, p2p_m


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
        print('wrong observatory!')
        sys.exit()


def get_winter_months():
    return np.array([1, 2, 3, 4, 12])


def get_summer_months():
    return np.array([7, 8, 9])


def get_all_months():
    return np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])


def get_intermediate_months():
    return np.array([5, 6, 10, 11])


def get_epoch(epoch):
    valid_epochs = {'winter': get_winter_months(),
                    'summer': get_summer_months(),
                    'all': get_all_months(),
                    'intermediate': get_intermediate_months()}
    return valid_epochs[epoch]


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


def get_gribfile_variables(file_name):
    """
    Function that returns all the different variable names in a grib file
    :param file_name:
    :return: varname (list): variable names
             varshortname (list): variable short names
    """
    print('opening file...')
    grb = pg.open(file_name)
    #grb.rewind()
    varshortname = []
    varname = []
    grb.read(10)[0]
    while True:
        v = grb.read(1)[0]
        vname = v.name.replace(" ", "")
        vshortname = v.shortName
        if vname not in varname:
            varname.append(vname)
            varshortname.append(vshortname)
        else:
            break
    grb.close()
    return varname, varshortname


def get_plevels(variable):
    plevels = []
    index = []
    for i in range(len(variable)):
        l = variable[i].level
        if l not in plevels:
            plevels.append(l)
            index.append(i+1)
        else:
            break
    return plevels, index


def get_grib_file_data(file_name):
    """
    This function opens a grib file, selects the parameters (all available: Temperature, Geopotential, RH, ...),
    and creates a dictionary variable where these parameters. It also returns the variable name (vn) and variable
    short name (vsn)

    Input: file_name (string)
           observatory (string). Possible values are 'north' or 'south'
           gridstep (float): grid spacing in degrees. Values are 1.0 for GDAS data and 0.75 for ECMWF data.

    Output: a txt file with the exact name as the input file name, but with .txt as extension
    """

    print('Working on %s' %(file_name))
    print('getting all variable names...')
    vn, vsn = get_gribfile_variables(file_name)
    print('indexing the file %s (this might take a while...)' % (file_name))
    grb = pg.index(file_name, 'shortName', 'typeOfLevel')
    print('selecting the parameters information for %s ...' % (file_name))
    data = []

    for sn in vsn:
        var = grb.select(shortName=sn, typeOfLevel='isobaricInhPa')
        data.append(var)
        gc.collect()

    datadict = dict(zip(vn, data))
    data = None
    gc.collect()
    return vn, vsn, datadict

def fill_RH_gaps(rhdata):
    """ in case the RH data for Plevels 20 and 50 hPa is not present, this function fills these gaps
        with 0.0"""
    RH = []
    for i in np.arange(len(rhdata)):
        RH.append(rhdata[i].values)
    RH = np.asarray(RH)

    for i in np.arange(len(RH)):
        if i*23 < len(RH):
            RH = np.insert(RH, (i*23, i*23), 0.0)
    return RH

def computedensity(p,T):
    return Ns * p / ps * Ts / T

def readgribfile2text(file_name, observatory, gridstep):
    """
    This function creates a txt file where the information from the get_grib_file_data function is written,
    together with date, year, month, day, hour, pressure level, real height and density.

    Input: file_name (string)
           observatory (string). Possible values are 'north' or 'south'
           gridstep (float): grid spacing in degrees. Values are 1.0 for GDAS data and 0.75 for ECMWF data.

    Output: a txt file with the exact name as the input file name, but with .txt as extension
    """

    if os.path.exists((file_name).split('.')[0] + '.txt'):
        print('Output file %s already exists. Aborting.' % ((file_name).split('.')[0] + '.txt'))
        sys.exit()

    vn, vsn, datadict = get_grib_file_data(file_name)

    latitude_obs, longitude_obs = get_observatory_coordinates(observatory)
    lat_gridpoint, lon_gridpoint = get_closest_gridpoint(latitude_obs, longitude_obs, gridstep)

    if len(datadict['Temperature']) != len(datadict['Relativehumidity']):
        RH = fill_RH_gaps(datadict['Relativehumidity'])
    else:
        RH = []
        for i in np.arange(len(datadict['Relativehumidity'])):
            RH.append(datadict['Relativehumidity'][i].values)
        RH = np.asarray(RH)

    # We create the table file and fill it with the information stored in the above variables, plus the height
    # and density computed form them.

    print('creating the txt file containing the selected data...')

    table_file = open(file_name.split('.')[0] + '.txt', 'w')
    print('Date year month day hour MJD Plevel T_average h n U V RH', file=table_file)

    for j in np.arange(len(datadict['Temperature'])):
        if (type(datadict['Temperature'][j].values) == float) or (len(datadict['Temperature'][j].values) == 1):
            if 'GeopotentialHeight' in vn:
                h = GetAltitudeFromGeopotentialHeight(datadict['GeopotentialHeight'][j].values, observatory)
            else:
                h = GetAltitudeFromGeopotential(datadict['Geopotential'][j].values, observatory)
            density = computedensity(datadict['Temperature'][j].level, datadict['Temperature'][j].values)
            mjd = date2mjd(datadict['Temperature'][j].year, datadict['Temperature'][j].month,
                           datadict['Temperature'][j].day, datadict['Temperature'][j].hour)
            print(int(datadict['Temperature'][j].dataDate), datadict['Temperature'][j].year,
                  datadict['Temperature'][j].month, datadict['Temperature'][j].day, datadict['Temperature'][j].hour, mjd,
                  datadict['Temperature'][j].level, datadict['Temperature'][j].values, h, density,
                  datadict['Ucomponentofwind'][j].values, datadict['Vcomponentofwind'][j].values, RH[j], file=table_file)

        else: # this is just in case the grib file contains more than one grid point
            if 'GeopotentialHeight' in vn:
                h = GetAltitudeFromGeopotentialHeight(np.float(datadict['GeopotentialHeight'][j].values[
                                 (datadict['GeopotentialHeight'][j].data()[1] == lat_gridpoint) &
                                 (datadict['GeopotentialHeight'][j].data()[2] == lon_gridpoint)]), observatory)
            else:
                h = GetAltitudeFromGeopotential(np.float(datadict['Geopotential'][j].values[
                                 (datadict['Geopotential'][j].data()[1] == lat_gridpoint) &
                                 (datadict['Geopotential'][j].data()[2] == lon_gridpoint)]), observatory)
            temperature = np.float(datadict['Temperature'][j].values[
                                       (datadict['Temperature'][j].data()[1] == lat_gridpoint) &
                                       (datadict['Temperature'][j].data()[2] == lon_gridpoint)])

            density = computedensity(datadict['Temperature'][j].level, temperature)
            mjd = date2mjd(datadict['Temperature'][j].year, datadict['Temperature'][j].month,
                           datadict['Temperature'][j].day, datadict['Temperature'][j].hour)

            print(int(datadict['Temperature'][j].dataDate), datadict['Temperature'][j].year, datadict['Temperature'][j].month,
                  datadict['Temperature'][j].day, datadict['Temperature'][j].hour, mjd, datadict['Temperature'][j].level,
                  temperature, h, density, datadict['Ucomponentofwind'][j].values,
                  datadict['Vcomponentofwind'][j].values, RH[j], file=table_file)

    table_file.close()
    datadict = None


def readgribfile2magic(file_name, observatory, gridstep):
    """
    This function opens a grib file, selects all parameters
    and finally creates a txt file where these parameters, together with date, year, month,
    day, hour, pressure level, real height and density, are written.
    Input: file_name (string)
           observatory (string). Possible values are 'north' or 'south'
           gridstep (float): grid spacing (0.75 degrees for ECMWF and 1.0 degrees for GDAS)
    Output: a txt file with the exact name as the input file name, but with .txt as extension in a format that can be
            read by MARS
    """
    if os.path.exists((file_name).split('.')[0] + 'MAGIC_format.txt'):
        print('Output file %s already exists. Aborting.' % ((file_name).split('.')[0] + 'MAGIC_format.txt'))
        sys.exit()

    vn, vsn, datadict = get_grib_file_data(file_name)

    pl, pl_index = get_plevels(datadict['Temperature'])
    new_pl_index = pl_index * int((len(datadict['Temperature'])/len(pl_index)))
    latitude_obs, longitude_obs = get_observatory_coordinates(observatory)
    lat_gridpoint, lon_gridpoint = get_closest_gridpoint(latitude_obs, longitude_obs, gridstep)

    # We create the table file and fill it with the information stored in the above variables, plus the height
    # and density computed form them.

    print('creating the txt file containing the selected data...')
    table_file = open(file_name.split('.')[0] + 'MAGIC_format.txt', 'w')

    for j in np.arange(len(datadict['Temperature'])):
        if (type(datadict['Temperature'][j].values) == float) or (len(datadict['Temperature'][j].values) == 1):
            if new_pl_index[j] == 1:
                print(str([0.00] * 34)[1:-1].replace(",", " "), file=table_file)
            if 'GeopotentialHeight' in vn:
                h = GetAltitudeFromGeopotentialHeight(datadict['GeopotentialHeight'][j].values, observatory)
            else:
                h = GetAltitudeFromGeopotential(datadict['Geopotential'][j].values, observatory)

            fields = (datadict['Temperature'][j].year - 2000, datadict['Temperature'][j].month,
                  datadict['Temperature'][j].day, datadict['Temperature'][j].hour, new_pl_index[j], h,
                  datadict['Temperature'][j].values, datadict['Ucomponentofwind'][j].values,
                  datadict['Vcomponentofwind'][j].values, RH[j])
            row_str = '{: >6d}{: >6d}{: >6d}{: >6d}{: >6d}{: >10.2f}{: >10.2f}{: >10.2f}{: >10.2f}{: >10.2f}'
            row_str = row_str.format(*fields)
            table_file.write(row_str + '\n')
        else:  # this is just in case the grib file contains more than one grid point
            if new_pl_index[j] == 1:
                print(str([0.00] * 34)[1:-1].replace(",", " "), file=table_file)
            if 'GeopotentialHight' in vn:
                h = np.float(datadict['GeopotentialHeight'][j].values[
                                 (datadict['GeopotentialHeight'][j].data()[1] == lat_gridpoint) &
                                 (datadict['GeopotentialHeight'][j].data()[2] == lon_gridpoint)])
            else:
                h = GetAltitudeFromGeopotential(datadict['Geopotential'][j].values, observatory)

            temperature = np.float(datadict['Temperature'][j].values[
                                       (datadict['Temperature'][j].data()[1] == lat_gridpoint) &
                                       (datadict['Temperature'][j].data()[2] == lon_gridpoint)])

            fields = (datadict['Temperature'][j].year - 2000, datadict['Temperature'][j].month,
                  datadict['Temperature'][j].day, datadict['Temperature'][j].hour, new_pl_index[j], h, temperature,
                  datadict['Ucomponentofwind'][j].values, datadict['Vcomponentofwind'][j].values,
                  RH[j])
            row_str = '{: >6d}{: >6d}{: >6d}{: >6d}{: >6d}{: >10.2f}{: >10.2f}{: >10.2f}{: >10.2f}{: >10.2f}'
            row_str = row_str.format(*fields)
            table_file.write(row_str + '\n')
    table_file.close()


def readgribfile2magic_fromtxt(txt_file):
#    TODO create the function
    """
    :param txt_file:
    :return:
    """
    return None


def runInParallel(function_name, list_of_gribfiles, observatory, gridstep):
    if multiprocessing.cpu_count() == 4:
        max_cpus = 2
    elif multiprocessing.cpu_count() >= 10:
        max_cpus = 5
    elif multiprocessing.cpu_count() == 1:
        max_cpus = 1
    else:
        max_cpus = multiprocessing.cpu_count() - 1

    first_element = 0
    while first_element + max_cpus <= len(list_of_gribfiles):
        sub_list_of_gribfiles = list_of_gribfiles[first_element:first_element + max_cpus]
        proc = []
        for f in sub_list_of_gribfiles:
            p = Process(target=function_name, args=(f, observatory, gridstep))
            proc.append(p)
            p.start()
        for p in proc:
            p.join()
        first_element += max_cpus
        if first_element + max_cpus > len(list_of_gribfiles):
            sub_list_of_gribfiles = list_of_gribfiles[first_element:]
            for f in sub_list_of_gribfiles:
                p = Process(target=readgribfile2text, args=(f, observatory, gridstep))
                proc.append(p)
                p.start()
            for p in proc:
                p.join()

def print_help():
    print("Usage: python grib_utils.py <options>")
    print("Options are:")
    print("        -r         <grib_file_name> <observatory> <gridstep>")
    print("                   note that <gridstep> is 0.75deg for ECMWF data")
    print("                   and 1.0 deg for GDAS data")
    print("        -rmagic    <grib_file_name> <observatory> <gridstep>")
    print("                   note that <gridstep> is 0.75deg for ECMWF data")
    print("                   and 1.0 deg for GDAS data")
    print("        -mjd       <mjd>")
    print("        -date      <yyyy-mm-dd-hh>")
    print(" ")
    print("                   Note: with the -r or -rmagic option, if a txt file")
    print("                   containing a list of grib files is passed instead")
    print("                   of a single grib file, the processing is run in parallel")
    print("                   using a certain number of CPU's")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == '-h' or sys.argv[1] == '-help':
        print_help()
        sys.exit()

    else:
        if sys.argv[1] == '-r':
            if sys.argv[2].split('.')[1] == 'grib' or sys.argv[2].split('.')[1] == 'grb':
                readgribfile2text(sys.argv[2], sys.argv[3], float(sys.argv[4]))
            elif sys.argv[2].split('.')[1] == 'txt' or sys.argv[2].split('.')[1] == 'dat':
                list_file = open(sys.argv[2], 'r')
                line = list_file.readline()
                list_of_files = []
                while line:
                    list_of_files.append(line[:-1])
                    line = list_file.readline()
                runInParallel(readgribfile2text, list_of_files, sys.argv[3], float(sys.argv[4]))

        elif sys.argv[1] == '-rmagic':
            if sys.argv[2].split('.')[1] == 'grib' or sys.argv[2].split('.')[1] == 'grb':
                readgribfile2magic(sys.argv[2], sys.argv[3], float(sys.argv[4]))
            elif sys.argv[2].split('.')[1] == 'txt' or sys.argv[2].split('.')[1] == 'dat':
                list_file = open(sys.argv[2], 'r')
                line = list_file.readline()
                list_of_files = []
                while line:
                    list_of_files.append(line[:-1])
                    line = list_file.readline()
                runInParallel(readgribfile2magic, list_of_files, sys.argv[3], float(sys.argv[4]))

        elif sys.argv[1] == '-mjd':
            print(mjd2date(float(sys.argv[2])))

        elif sys.argv[1] == '-date':
            date = sys.argv[2].split('-')
            print(date2mjd(int(date[0]), int(date[1]), int(date[2]), int(date[3])))

        else:
            print('Wrong option...')
            print_help()
            sys.exit()

