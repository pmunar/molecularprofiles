from tqdm import *
from molecularprofiles.utils.grib_utils import date2mjd, get_epoch
import numpy as np


def read_file(file, epoch_text):
    global months, month
    print("loading and selecting data")
    file = open(file)
    date, year, month, day, hour, p, T, h, n, U, V, RH = np.loadtxt(file, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                                                                   11), unpack=True)

    mjd = []

    for i in tqdm(np.arange(len(date))):
        mjd.append(date2mjd(year[i], month[i], day[i], hour[i]))
    print('\n')
    mjd = np.asarray(mjd)
    epoch = get_epoch(epoch_text)

    # mjd = mjd[month == epoch]
    # year = year[month == epoch]
    # day = date[month == epoch]
    # hour = hour[month == epoch]
    # h   = h[month == epoch]
    # dn   = n[month == epoch]
    # p   = p[month == epoch]
    # month = month[month == epoch]

    if epoch_text == 'winter':
        mjd = mjd[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                              (month == epoch[3]) | (month == epoch[4])]
        h = h[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3]) | (month == epoch[4])]
        n = n[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3]) | (month == epoch[4])]
        p = p[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3]) | (month == epoch[4])]
        T = T[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3]) | (month == epoch[4])]
        U = U[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3]) | (month == epoch[4])]
        V = V[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3]) | (month == epoch[4])]
        RH = RH[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3]) | (month == epoch[4])]

    elif epoch_text == 'summer':
        mjd = mjd[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        h = h[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        n = n[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        p = p[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        T = T[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        U = U[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        V = V[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        RH = RH[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]

    elif epoch_text == 'intermediate':
        mjd = mjd[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                              (month == epoch[3])]
        h = h[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3])]
        n = n[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3])]
        p = p[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                          (month == epoch[3])]
        T = T[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
              (month == epoch[3])]
        U = U[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
              (month == epoch[3])]
        V = V[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
              (month == epoch[3])]
        RH = RH[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
              (month == epoch[3])]

    elif epoch_text == 'all':
        mjd = mjd
        h = h
        n = n
        p = p
        T = T
        U = U
        V = V
        RH = RH



    return mjd, year, month, day, hour, p, h, n, T, U, V, RH
