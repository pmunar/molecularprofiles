from tqdm import *
from molecularprofiles.utils.grib_utils import date2mjd, get_epoch
import numpy as np


def read_ecmwf(file_ecmwf, epoch_text):
    global months, month_ecmwf
    print("loading and selecting ecmwf data")
    ecmwf_file = open(file_ecmwf)
    date_ecmwf, year_ecmwf, month_ecmwf, day_ecmwf, hour_ecmwf, p_ecmwf, T_ecmwf, h_ecmwf, n_ecmwf, U_ecmwf, V_ecmwf, \
        RH_ecmwf = np.loadtxt(ecmwf_file, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True)

    mjd_ecmwf = []

    for i in tqdm(np.arange(len(date_ecmwf))):
        mjd_ecmwf.append(date2mjd(year_ecmwf[i], month_ecmwf[i], day_ecmwf[i], hour_ecmwf[i]))
    print('\n')
    mjd_ecmwf = np.asarray(mjd_ecmwf)
    epoch = get_epoch(epoch_text)

    # mjd_ecmwf = mjd_ecmwf[month_ecmwf == epoch]
    # year_ecmwf = year_ecmwf[month_ecmwf == epoch]
    # day_ecmwf = date_ecmwf[month_ecmwf == epoch]
    # hour_ecmwf = hour_ecmwf[month_ecmwf == epoch]
    # h_ecmwf   = h_ecmwf[month_ecmwf == epoch]
    # dn_ecmwf   = n_ecmwf[month_ecmwf == epoch]
    # p_ecmwf   = p_ecmwf[month_ecmwf == epoch]
    # month_ecmwf = month_ecmwf[month_ecmwf == epoch]

    if epoch_text == 'winter':
        mjd_ecmwf = mjd_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                              (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]
        h_ecmwf = h_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]
        n_ecmwf = n_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]
        p_ecmwf = p_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]

    elif epoch_text == 'summer':
        mjd_ecmwf = mjd_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]
        h_ecmwf = h_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]
        n_ecmwf = n_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]
        p_ecmwf = p_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]

    elif epoch_text == 'intermediate':
        mjd_ecmwf = mjd_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                              (month_ecmwf == epoch[3])]
        h_ecmwf = h_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3])]
        n_ecmwf = n_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3])]
        p_ecmwf = p_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3])]

    elif epoch_text == 'all':
        mjd_ecmwf = mjd_ecmwf
        h_ecmwf = h_ecmwf
        n_ecmwf = n_ecmwf
        p_ecmwf = p_ecmwf



    return mjd_ecmwf, year_ecmwf, month_ecmwf, day_ecmwf, hour_ecmwf, h_ecmwf, p_ecmwf, n_ecmwf