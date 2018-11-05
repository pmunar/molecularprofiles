from tqdm import *
from molecularprofiles.utils.grib_utils import date2mjd, get_epoch
import numpy as np


def select_dataframe_epoch(file, epoch_text):
    global months, month
    print("loading and selecting data")
    file = open(file)
    date, year, month, day, hour, mjd, p, T, h, n, U, V, RH = np.loadtxt(file, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                                        10,11, 12), unpack=True)

    epoch = get_epoch(epoch_text)

    if epoch_text == 'winter':
        date = date[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                  (month == epoch[3]) | (month == epoch[4])]
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
        date = date[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        mjd = mjd[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        h = h[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        n = n[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        p = p[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        T = T[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        U = U[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        V = V[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]
        RH = RH[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2])]

    elif epoch_text == 'intermediate':
        date = date[(month == epoch[0]) | (month == epoch[1]) | (month == epoch[2]) |
                  (month == epoch[3])]
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
