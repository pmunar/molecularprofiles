#!/usr/bin/env python

from builtins import str
import numpy as np
import pygrib as pg
from molecularprofiles.utils.meteorological_constants import *
from molecularprofiles.utils.observatory import *
import sys
import gc
import multiprocessing
from multiprocessing import Process
import os
from tqdm import tqdm
import argparse



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


def GetAltitudeFromGeopotential(geop_height, latitude_obs):
    """
    Function to compute the real altitude from the geopotential value at a certain coordinates on Earth
    :param geop_height:
    :param latitude_obs: geographical latitude of interest in degrees
    :return: real altitude as fGeoidOffset (in m)
    """
    latitude = np.radians(latitude_obs)
    geop_heightkm = geop_height / 1000. / 9.80665  # dividing by the acceleration of gravity on Earth
    cos2lat = np.cos(2 * latitude)
    # convert from geopotential height to geometric altitude:
    # Uses expression 20 from http://www.earth.sinica.edu.tw/~bfchao/publication/eng/2005-Boy&Chao_JGR2005_precise%20evaluation%20of%20atmospheric%20loading%20effects%20on%20earths%20time-variable%20gravity%20field.pdf
    z = (1. + 0.002644 * cos2lat) * geop_heightkm + (1 + 0.0089 * cos2lat) * geop_heightkm * geop_heightkm / 6245.
    # convert Z to meter
    return 1.0E3 * z  # This is fGeoidOffset


def GetAltitudeFromGeopotentialHeight(geop, latitude_obs):
    """
    Function to compute the real altitude from the geopotential value at a certain coordinates on Earth
    :param geop_height:
    :param latitude_obs: geographical latitude of interest in degrees
    :return: real altitude as fGeoidOffset (in m)
    """
    latitude = np.radians(latitude_obs)
    geop_km = geop / 1000.
    cos2lat = np.cos(2 * latitude)
    # convert from geopotential height to geometric altitude:
    # Uses expression 20 from http://www.earth.sinica.edu.tw/~bfchao/publication/eng/2005-Boy&Chao_JGR2005_precise%20evaluation%20of%20atmospheric%20loading%20effects%20on%20earths%20time-variable%20gravity%20field.pdf
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
    and creates a dictionary variable with these parameters. It also returns the variable name (vn) and variable
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
    print('selecting the parameters information for %s (this might take a while...)' % (file_name))
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


def compute_wind_direction(u,v):
    angle = np.arctan2(-1*u,-1*v)*180./np.pi
    if angle < 0.:
        angle += 360.
    direction = angle
    return direction


def compute_wind_speed(u,v):
    return np.sqrt(u**2. + v**2.)


def readgribfile2text(file_name, gridstep, observatory=None, lat=None, lon=None):
    """
    This function creates a txt file where the information from the get_grib_file_data function is written,
    together with date, year, month, day, hour, pressure level, real height and density.

    Input: file_name (string)
           observatory (string). Possible values are 'north', 'south' or any other name. If NOT north or south,
           then the program asks for the coordinates.
           gridstep (float): grid spacing in degrees. Values are 1.0 for GDAS data and 0.75 for ECMWF data.
           lat: (float, optional) latitude of the observatory in degrees
           lon: (float, optional) longitude of the observatory in degrees
    Output: a txt file with the exact name as the input file name, but with .txt as extension
    """

    if os.path.exists(os.path.splitext(file_name)[0] + '.txt'):
        print('Output file %s already exists. Aborting.' % (os.path.splitext(file_name)[0] + '.txt'))
        sys.exit()

    vn, vsn, datadict = get_grib_file_data(file_name)

    latitude_obs, longitude_obs = None, None

    if observatory:
        latitude_obs, longitude_obs = get_observatory_coordinates(observatory)
    elif lat and lon:
        latitude_obs, longitude_obs = lat, lon

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
    print('Date year month day hour MJD P Temp h n n_exp U V wind_speed wind_direction RH', file=table_file)

    pbar = tqdm(total=len(datadict['Temperature']))
    for j in np.arange(len(datadict['Temperature'])):
        pbar.update(1)
        if (type(datadict['Temperature'][j].values) == float) or (len(datadict['Temperature'][j].values) == 1):
            if 'GeopotentialHeight' in vn:
                h = GetAltitudeFromGeopotentialHeight(datadict['GeopotentialHeight'][j].values, latitude_obs)
            else:
                h = GetAltitudeFromGeopotential(datadict['Geopotential'][j].values, latitude_obs)
            density = computedensity(datadict['Temperature'][j].level, datadict['Temperature'][j].values)
            density_exp = density / Ns * np.exp(h/Hs)
            wind_speed = compute_wind_speed(datadict['Ucomponentofwind'][j].values, datadict['Vcomponentofwind'][j].values)
            wind_direction = compute_wind_direction(datadict['Ucomponentofwind'][j].values, datadict['Vcomponentofwind'][j].values)
            mjd = date2mjd(datadict['Temperature'][j].year, datadict['Temperature'][j].month,
                           datadict['Temperature'][j].day, datadict['Temperature'][j].hour)
            print(int(datadict['Temperature'][j].dataDate), datadict['Temperature'][j].year,
                  datadict['Temperature'][j].month, datadict['Temperature'][j].day, datadict['Temperature'][j].hour, mjd,
                  datadict['Temperature'][j].level, datadict['Temperature'][j].values, h, density, density_exp,
                  datadict['Ucomponentofwind'][j].values, datadict['Vcomponentofwind'][j].values, wind_speed,
                  wind_direction, RH[j], file=table_file)

        else: # this is just in case the grib file contains more than one grid point
            if 'GeopotentialHeight' in vn:
                h = GetAltitudeFromGeopotentialHeight(np.float(datadict['GeopotentialHeight'][j].values[
                                 (datadict['GeopotentialHeight'][j].data()[1] == lat_gridpoint) &
                                 (datadict['GeopotentialHeight'][j].data()[2] == lon_gridpoint)]), latitude_obs)
            else:
                h = GetAltitudeFromGeopotential(np.float(datadict['Geopotential'][j].values[
                                 (datadict['Geopotential'][j].data()[1] == lat_gridpoint) &
                                 (datadict['Geopotential'][j].data()[2] == lon_gridpoint)]), latitude_obs)
            temperature = np.float(datadict['Temperature'][j].values[
                                       (datadict['Temperature'][j].data()[1] == lat_gridpoint) &
                                       (datadict['Temperature'][j].data()[2] == lon_gridpoint)])

            density = computedensity(datadict['Temperature'][j].level, temperature)
            wind_speed = compute_wind_speed(datadict['Ucomponentofwind'][j].values, datadict['Vcomponentofwind'][j].values)
            wind_direction = compute_wind_direction(datadict['Ucomponentofwind'][j].values, datadict['Vcomponentofwind'][j].values)
            density_exp = density/Ns * np.exp(h/Hs)
            mjd = date2mjd(datadict['Temperature'][j].year, datadict['Temperature'][j].month,
                           datadict['Temperature'][j].day, datadict['Temperature'][j].hour)

            print(int(datadict['Temperature'][j].dataDate), datadict['Temperature'][j].year, datadict['Temperature'][j].month,
                  datadict['Temperature'][j].day, datadict['Temperature'][j].hour, mjd, datadict['Temperature'][j].level,
                  temperature, h, density, density_exp, datadict['Ucomponentofwind'][j].values,
                  datadict['Vcomponentofwind'][j].values, wind_speed, wind_direction, RH[j], file=table_file)

    table_file.close()
    pbar.close()
    datadict = None


def readgribfile2magic(file_name, gridstep, observatory=None, lat=None, lon=None):
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
    if os.path.exists(os.path.splitext(file_name)[0] + '.txt'):
        print('Output file %s already exists. Aborting.' % (os.path.splitext(file_name)[0] + '.txt'))
        sys.exit()

    vn, vsn, datadict = get_grib_file_data(file_name)

    pl, pl_index = get_plevels(datadict['Temperature'])
    new_pl_index = pl_index[::-1] * int((len(datadict['Temperature'])/len(pl_index)))

    if observatory:
        latitude_obs, longitude_obs = get_observatory_coordinates(observatory)
    else:
        latitude_obs, longitude_obs = lat, lon

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
                h = GetAltitudeFromGeopotentialHeight(datadict['GeopotentialHeight'][j].values, latitude_obs)
            else:
                h = GetAltitudeFromGeopotential(datadict['Geopotential'][j].values, latitude_obs)

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
                h = GetAltitudeFromGeopotential(datadict['Geopotential'][j].values, latitude_obs)

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
    """
    :param txt_file:
    :return:
    """
    input_f = open(txt_file, 'r')
    output_f = open(txt_file + 'MAGIC_format.txt', 'w')

    date, year, month, day, hour, mjd, p, T, h, n, U, V, RH = np.loadtxt(input_f, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12)
                                                                         , unpack=True, skiprows=1)

    pl = np.unique(p)[::-1]
    pl_index = (np.arange(len(pl))+1).tolist()
    new_pl_index = pl_index * int((len(T) / len(pl_index)))
    remaining_index = len(p) - len(new_pl_index)
    new_pl_index = new_pl_index + pl_index[: remaining_index]

    pbar = tqdm(total=len(T))
    for j in np.arange(len(T)):
        pbar.update(1)
        if new_pl_index[j] == 1:
            print(str([0.00] * 34)[1:-1].replace(",", " "), file=output_f)

        fields = (int(year[j] - 2000), int(month[j]), int(day[j]), int(hour[j]), int(new_pl_index[j]), h[j], T[j], U[j],
                  V[j], RH[j])
        row_str = '{: >6d}{: >6d}{: >6d}{: >6d}{: >6d}{: >10.2f}{: >10.2f}{: >10.2f}{: >10.2f}{: >10.2f}'
        row_str = row_str.format(*fields)
        output_f.write(row_str + '\n')

    pbar.close()

def runInParallel(function_name, file_list, gridstep, observatory=None, lat=None, lon=None):
    if multiprocessing.cpu_count() == 4:
        max_cpus = 2
    elif multiprocessing.cpu_count() >= 10:
        max_cpus = 10
    elif multiprocessing.cpu_count() == 1:
        max_cpus = 1
    else:
        max_cpus = multiprocessing.cpu_count() - 1
    fname = open(file_list)
    line = fname.readline()
    list_of_gribfiles = []
    while line:
        list_of_gribfiles.append(line[:-1])
        line = fname.readline()

    first_element = 0
    while first_element + max_cpus <= len(list_of_gribfiles):
        sub_list_of_gribfiles = list_of_gribfiles[first_element:first_element + max_cpus]
        proc = []
        for f in sub_list_of_gribfiles:
            if observatory:
                p = Process(target=function_name, args=(f, gridstep, observatory, lat, lon))
            elif lat and lon:
                p = Process(target=function_name, args=(f, gridstep, observatory, lat, lon))
            proc.append(p)
            p.start()
        for p in proc:
            p.join()
        first_element += max_cpus
        if first_element + max_cpus > len(list_of_gribfiles):
            sub_list_of_gribfiles = list_of_gribfiles[first_element:]
            for f in sub_list_of_gribfiles:
                if observatory:
                    p = Process(target=function_name, args=(f, gridstep, observatory, lat, lon))
                elif lat and lon:
                    p = Process(target=function_name, args=(f, gridstep, observatory, lat, lon))
                proc.append(p)
                p.start()
            for p in proc:
                p.join()

def merge_txt_from_grib(txtfile, output_file='merged_from_grib.txt'):
    lf = open(txtfile, 'r')
    outfile = open(output_file, 'w')

    line = lf.readline()
    first = True
    while line:
        datafile = open(line[:-1], 'r')
        if first:
            dataline = datafile.readline()
        else:
            datafile.readline()
            dataline = datafile.readline()

        while dataline:
            print(dataline[:-1], file=outfile)
            dataline = datafile.readline()
        first = False
        datafile.close()
        line = lf.readline()
    lf.close()
    outfile.close()


def print_help():
    print("Usage: python grib_utils.py <options>")
    print("Options are:")
    print("        -r         <grib_file_name> <gridstep> <observatory>")
    print("                   note that <gridstep> is 0.75deg for ECMWF data")
    print("                   and 1.0 deg for GDAS data")
    print("        -rmagic    <grib_file_name> <observatory> <gridstep>")
    print("                   note that <gridstep> is 0.75deg for ECMWF data")
    print("                   and 1.0 deg for GDAS data")
    print("        -mjd       <mjd>")
    print("        -date      <yyyy-mm-dd-hh>")
    print("        -merge or -m <list_of_txt_files> <output_name>")
    print(" ")
    print("                   Note: with the -r or -rmagic option, if a txt file")
    print("                   containing a list of grib files is passed instead")
    print("                   of a single grib file, the processing is run in parallel")
    print("                   using a certain number of CPU's")


#if __name__ == "__main__":

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--grib_file', help='the grib file to process')
parser.add_argument('-gridstep', help='the gridstep in degrees. If GDAS or GFS data, gridstep=1.0 deg; '
                                                   'If ECMWF data, gridstep=0.75 deg', type=float)
parser.add_argument('-o', '--observatory', help='north or south. If no observatory is provided, '
                                                                 'then the system asks for the coordinates of interest')
parser.add_argument('-c', '--coordinates', nargs=2, help='latitude and longitude of the place '
                                                                               'of interest, in degrees', type=float)
parser.add_argument('-r', action='store_true', help='<grib_file_name> <gridstep> <observatory> \n '
                               'note that <gridstep> is 0.75deg for ECMWF data \n'
                               'and 1.0 deg for GDAS data. If a txt file containing \n '
                               'a list of grib files is passed instead of a single \n '
                               'grib file, the processing is run in parallel \n'
                               'using a certain number of CPUs')
parser.add_argument('-rmagic', action='store_true', help='<grib_file_name> <observatory> <gridstep> \n '
                               'note that <gridstep> is 0.75deg for ECMWF data \n '
                               'and 1.0 deg for GDAS data.  If a txt file containing \n '
                               'a list of grib files is passed instead of a single \n '
                               'grib file, the processing is run in parallel \n'
                               'using a certain number of CPUs')
parser.add_argument('--parallel', action='store_true', help='if a txt file containing a list of grib files is passed as -g option,'
                                             'it opens them and extracts information in parallel, one grib file per'
                                             'CPU')
parser.add_argument('-mjd', help= 'if selected, transforms MJD information into date in YYYY MM DD HH format', type=float)
parser.add_argument('-date', help='if selected, transforms date information into MJD format', type=str)
parser.add_argument('-m', '--merge', nargs='+', help='followed by a filename containing a list of txt files, \n '
                                   'it merges them into a single txt file')


if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    if args.grib_file:
        if args.observatory:
            if args.grib_file.lower().endswith(('.grib','.grb', '.txt', '.dat')):
                if args.r and not args.parallel:
                    readgribfile2text(args.grib_file, args.gridstep, observatory=args.observatory)
                elif args.rmagic:
                    readgribfile2magic(args.grib_file, args.gridstep, observatory=args.observatory)
                elif args.r and args.parallel:
                    print('EXECUTING THIS')
                    runInParallel(readgribfile2text, args.grib_file, args.gridstep, observatory=args.observatory)
                elif args.rmagic and args.parallel:
                    runInParallel(readgribfile2magic, args.grib_file, args.gridstep, observatory=args.observatory)
            else:
                print('file extension not recognized. Exiting')
                sys.exit()
        elif not args.observatory and args.coordinates:
            lat, lon = args.coordinates
            if args.grib_file.lower().endswith(('.grib','.grb', '.txt', '.dat')):
                if args.r:
                    readgribfile2text(args.grib_file, args.gridstep, lat=lat, lon=lon)
                elif args.rmagic:
                    readgribfile2magic(args.grib_file, args.gridstep, lat=lat, lon=lon)
                elif args.r and args.parallel:
                    runInParallel(readgribfile2text, args.grib_file, args.gridstep, lat=lat, lon=lon)
                elif args.rmagic and args.parallel:
                    runInParallel(readgribfile2magic, args.grib_file, args.gridstep, lat=lat, lon=lon)
            else:
                print('file extension not recognized. Exiting')
                sys.exit()
        else:
            print('Too many options. Please specify either observatory only, or coordinates only')

    elif args.mjd:
        print(mjd2date(args.mjd))

    elif args.date:
        date = args.date.split('-')
        print(date2mjd(int(date[0]), int(date[1]), int(date[2]), int(date[3])))

    elif args.merge:
        if isinstance(args.merge, list) and len(args.merge) == 2:
            merge_txt_from_grib(args.merge[0], output_file=args.merge[1])
        else:
            merge_txt_from_grib(args.merge[0])