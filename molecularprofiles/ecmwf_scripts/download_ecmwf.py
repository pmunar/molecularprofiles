#!/usr/bin/env python

import numpy as np
from ecmwfapi import ECMWFDataServer
import datetime
import sys
import calendar
import argparse

def retrieve_interim(date_start, date_end, latitude, longitude, days=7, outtag='my_ecmwf_interim_data'):
    """
       A function to demonstrate how to iterate efficiently over several years and months etc
       for a particular interim_request.
       Change the variables below to adapt the iteration to your needs.
       You can use the variable 'target' to organise the requested data in files as you wish.
       In the example below the data are organised in files per month. (eg "interim_daily_201510.grb")

        inputs:
            date_start : (string) starting date of the period to download in "YYYY-MM-DD" format
            date_end   : (string) ending date of the period to download in "YYYY-MM-DD" format
            latitude   : (float) geographical latitude of the place of interest in degrees
            longitude  : (float) geographical longitude of the place of interest in degrees.
                        If negative it will be transformed to positive.
            outtag (optional) : (string) tag to append to the output file name

        outputs:
            outfile : grib2 format files, as many as months from date_start until date_end
            (e.g. if date_start = "2017-01-01" and date_end = "2017-03-05" it will produce
            three files, one containing all January data, one containing all February data
            and another one containing the first 5 days of March data)

    """
    date_start = list(map(int, date_start.split('-')))
    date_end = list(map(int, date_end.split('-')))
    yearStart, monthStart, dayStart = date_start
    yearEnd, monthEnd, dayEnd = date_end

    date_start = datetime.datetime(yearStart, monthStart, dayStart)
    date_end = datetime.datetime(yearEnd, monthEnd, dayEnd)

    if dayStart > calendar.monthrange(yearStart, monthStart)[1] or \
                    dayEnd > calendar.monthrange(yearEnd, monthEnd)[1]:
        print('Wrong date! Check input dates')
        exit()
    if date_start > date_end:
        print('Wrong date. Date Start must be prior to Date End')

    while True:
        date_end_e = date_start + datetime.timedelta(days=days-1)
        if date_end_e < date_end:
            start_date = '%04d%02d%02d' % (date_start.year, date_start.month, date_start.day)
            end_date = '%04d%02d%02d' % (date_end_e.year, date_end_e.month, date_end_e.day)
            print('Downloading data between {} and {}'.format(start_date, end_date))
            outfile = '{}_{}_to_{}.grib'.format(outtag, start_date, end_date)
            request_ecwmf(start_date, end_date, latitude, longitude, outfile)
            date_start = date_end_e
        else:
            date_end_e = date_end
            start_date = '%04d%02d%02d' % (date_start.year, date_start.month, date_start.day)
            end_date = '%04d%02d%02d' % (date_end_e.year, date_end_e.month, date_end_e.day)
            print('Downloading data between {} and {}'.format(start_date, end_date))
            outfile = '{}_{}_to_{}.grib'.format(outtag, start_date, end_date)
            request_ecwmf(start_date, end_date, latitude, longitude, outfile)
            break

    # date_start_split = date_start.split('-')
    # yearStart = int(date_start_split[0])
    # monthStart = int(date_start_split[1])
    # dayStart = int(date_start_split[2])
    #
    # date_end_split = date_end.split('-')
    # yearEnd = int(date_end_split[0])
    # monthEnd = int(date_end_split[1])
    # dayEnd = int(date_end_split[2])
    #
    # if dayStart > calendar.monthrange(yearStart, monthStart)[1] or \
    #                 dayEnd > calendar.monthrange(yearEnd, monthEnd)[1]:
    #     print('Wrong date! Check input dates')
    #     exit()
    #
    # first = True
    # last = False
    # for year in list(range(yearStart, yearEnd + 1)):
    #     if first:
    #         month = monthStart
    #     else:
    #         month = 1
    #     while last == False and month < 13:
    #     #for month in list(range(monthStart, monthEnd + 1)):
    #         startDate = '%04d%02d%02d' % (year, month, dayStart)
    #         if year == yearEnd and month == monthEnd:
    #             last = True
    #             lastDate = '%04d%02d%02d' % (year, month, dayEnd)
    #             outfile = outtag + '_%04d%02d.grib' % (year, month)
    #             request_ecwmf(startDate, lastDate, latitude, longitude, outfile)
    #             break
    #         else:
    #             lastDate = '%04d%02d%02d' % (year, month, calendar.monthrange(year, month)[1])
    #         dayStart = 1
    #         outfile = outtag+'_%04d%02d.grib' % (year, month)
    #         request_ecwmf(startDate, lastDate, latitude, longitude, outfile)
    #         first = False
    #        month += 1


def request_ecwmf(date_i, date_f, lat, lon, outfile='my_ecmwf_file.grib'):

    """
    This function downloads ERA-Interim forecasts, on pressure levels.
    The data volume for all pressure level data is about 5GB per day,
    and all pressure level data for more than a single day will exceed
    the WebAPI limit of 600.000 fields. Thus please restrict the download
    to what you really need.
    It is likely you need to split your request, this is best done by time
    periods, ie. first download for month 1, then for month 2, etc.

    Inputs:
        date_i : string, initial date in string format as "YYYY-MM-DD"
        date_f : string, final date in string format as "YYYY-MM-DD"
        lat : float, latidude of interest in degrees
        lon : float, longitude of interest in degrees
        outfile : string, output file name

    Output:
        outfile : grib2 format file containing the downloaded data
    """

    server = ECMWFDataServer()

    step = 0.75  # grid step in degrees
    lons_grid = np.zeros(int(360 / step) + 1)
    lats_grid = np.zeros(int(180 / step) + 1)

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

    server.retrieve({
        # Specify the ERA-Interim data archive. Don't change.
        "class": "ei",
        "dataset": "interim",
        "expver": "1",
        "stream": "oper",
        # pressure levels (levtype:pl), all available levels (levelist)
        "levtype": "pl",
        "levelist": "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/\
                    800/825/850/875/900/925/950/975/1000",
        # forecast (type:fc), from both daily forecast runs (time) with all available forecast steps (step, in hours)
        "type": "an",
        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
        "step": "0",
        # all available parameters, for codes see http://apps.ecmwf.int/codes/grib/param-db
        "param": "60.128/129.128/130.128/131.128/132.128/133.128/135.128/138.128/155.128/157.128/203.128/246.128/\
                  247.128/248.128",
        # two days worth of data
        "date": date_i+'/to/'+date_f,
        # in 0.75 degrees lat/lon
        "grid": "0.75/0.75",
        # optionally restrict area to Europe (in N/W/S/E).
        "area": str(nearest_lat)+'/'+str(nearest_lon)+'/'+str(nearest_lat)+'/'+str(nearest_lon),
        # Optionally get output in NetCDF format. However, this does not work with the above request due to
        # overlapping timestamps.
        # "format" : "netcdf",
        # set an output file name
        "target": outfile,
    })


def find_nearest(a, num):  # Function to find the nearest grid position to a given latitude or longitude
    return a[np.abs(a-num).argmin()]


parser = argparse.ArgumentParser()
parser.add_argument('date start', help='the starting date of the dataset to download, in YYYY-MM-DD format')
parser.add_argument('date end', help='the ending date of the dataset to download, in YYYY-MM-DD format')
parser.add_argument('latitude', help='the latitude of the location of interest, in degrees')
parser.add_argument('longitude', help='the longitude of the location of interest, in degrees')

parser.add_argument('-d', '--days', help='number of days contained in each data file to download. '
                                         'Usual values are 7 or 30 days. Best value for RAM usage '
                                         'is 7 days')
parser.add_argument('-t', '--tag', help='the output tag to put in the downloaded data files')

if __name__ == '__main__':
    args = parser.parse_args()
    if not args.days and not args.tag:
        retrieve_interim(sys.argv[1], sys.argv[2], np.float(sys.argv[3]), np.float(sys.argv[4]))
    elif args.days:
        retrieve_interim(sys.argv[1], sys.argv[2], np.float(sys.argv[3]), np.float(sys.argv[4]), days=args.days)
    elif args.tag:
        retrieve_interim(sys.argv[1], sys.argv[2], np.float(sys.argv[3]), np.float(sys.argv[4]), outtag=args.tag)
    elif args.tag and args.days:
        retrieve_interim(sys.argv[1], sys.argv[2], np.float(sys.argv[3]), np.float(sys.argv[4]), outtag=args.tag,
                         days=args.days)
