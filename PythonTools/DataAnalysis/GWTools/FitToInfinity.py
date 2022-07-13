
"""
Use radii of extraction to make a fit and extrapolate to infinity

This uses a 1/r and a 1/r^2 fit to estimate the error
(see (3) and (4) of https://arxiv.org/pdf/1012.3173.pdf)

This is better than Richardson Extrapolation, because:
1. uses all the extraction radii, not just 2
2. is much much less sensitive to noise

"""

import sys, os, glob
sys.path.append("../Base/")
sys.path.append("../")
from utils import *
import SmallDataIOReader

import multiprocessing
from scipy.interpolate import make_interp_spline
import numpy as np

def data_fit_to_infinity_from_file(file_path, mass, dest = None):
    file = SmallDataIOReader.File(file_path)

    radii = get_radii(file_path, verbose=False)
    radii_indices = np.array((range(len(radii))))*2+1
    rs_extraction = [tortoise_radius(radius, mass) for radius in radii]

    time = file.getData()[:,0] / mass
    datas_re = [file.getData()[:,index] for index in radii_indices]
    datas_im = [file.getData()[:,index+1] for index in radii_indices]

    time_new = np.linspace(np.min(time)-rs_extraction[0]/mass, np.max(time)-rs_extraction[-1]/mass, 10000) 
    datas_re_spl = [make_interp_spline(time, data, k=3) for data in datas_re]
    datas_im_spl = [make_interp_spline(time, data, k=3) for data in datas_im]
    datas_re_smooth = np.array([data_spl(time_new + rs_extraction[i]/mass) for i, data_spl in enumerate(datas_re_spl)])
    datas_im_smooth = np.array([data_spl(time_new + rs_extraction[i]/mass) for i, data_spl in enumerate(datas_im_spl)])

    # FIT
    def fit_func_r2(p, rs):
        m, a, b = p
        return m + a/rs + b/rs**2
    def fit_func_r(p, rs):
        m, a = p
        return m + a/rs

    print(rs_extraction[0:])
    print(datas_re_smooth[0:,int(len(time_new)/2)])

    # fit to infinity
    start = 0
    data_re_fit_infinity_full_r2 = [FIT(rs_extraction[start:], datas_re_smooth[start:,time], fit_func_r2, 3,
                                    beta0=[datas_re_smooth[0,time], 0, 0], verbose=False) for time in range(len(time_new))]
    data_im_fit_infinity_full_r2 = [FIT(rs_extraction[start:], datas_im_smooth[start:,time], fit_func_r2, 3,
                                    beta0=[datas_im_smooth[0,time], 0, 0], verbose=False) for time in range(len(time_new))]
    data_re_fit_infinity_full_r = [FIT(rs_extraction[start:], datas_re_smooth[start:,time], fit_func_r, 2,
                                    beta0=[datas_re_smooth[0,time], 0], verbose=False) for time in range(len(time_new))]
    data_im_fit_infinity_full_r = [FIT(rs_extraction[start:], datas_im_smooth[start:,time], fit_func_r, 2,
                                    beta0=[datas_im_smooth[0,time], 0], verbose=False) for time in range(len(time_new))]

    print(data_re_fit_infinity_full_r2[int(len(time_new)/2)][0][0])
    print(data_re_fit_infinity_full_r[int(len(time_new)/2)][0][0])

    data_re_fit_infinity = [fit_r2[0][0] for fit_r2 in data_re_fit_infinity_full_r]
    data_im_fit_infinity = [fit_r2[0][0] for fit_r2 in data_im_fit_infinity_full_r]
    data_re_fit_errors = [max(fit_r2[1][0], fit_r[1][0], abs(fit_r2[0][0] - fit_r[0][0]))
                            for fit_r2, fit_r in zip(data_re_fit_infinity_full_r2, data_re_fit_infinity_full_r)]
    data_im_fit_errors = [max(fit_r2[1][0], fit_r[1][0], abs(fit_r2[0][0] - fit_r[0][0]))
                            for fit_r2, fit_r in zip(data_im_fit_infinity_full_r2, data_im_fit_infinity_full_r)]

    data_re_fit_first_radius_errors = [abs(fit_r2[0][0] - data)
                            for fit_r2, data in zip(data_re_fit_infinity_full_r2, datas_re_smooth[0])]
    data_im_fit_first_radius_errors = [abs(fit_r2[0][0] - data)
                            for fit_r2, data in zip(data_im_fit_infinity_full_r2, datas_im_smooth[0])]

    new_path = os.path.splitext(file_path if dest is None else dest)[0] + "_at_infinity.dat"
    np.savetxt(new_path,
        np.transpose([time_new, datas_re_smooth[0], datas_im_smooth[0], data_re_fit_infinity, data_im_fit_infinity, data_re_fit_first_radius_errors, data_im_fit_first_radius_errors,
            data_re_fit_errors, data_im_fit_errors]),
        header=f"(time-r_ex^*)/M \tRe_part_r{radii[0]} \tIm_part_r{radii[0]} \tRe_part_infty \tIm_part_infty \tRe_part_r{radii[0]}_err \tIm_part_r{radii[0]}_err \tRe_part_infty_err \tIm_part_infty_err")


################################################################################
# FIT TO INFINITY

# Below is just a fancy way of running the function above for many folders, in parallel

# pairs of folder + spacetime mass of the respective run
folders_for_infinity_fit = [("gm1", 1)]

data_folder = "/home/taigofr/Documents/MyDocuments/Education/PhD/Projects/EsGB_Llibert/" # folder where to look for folders

# filename in each folder
filename = "strain_??_cutoff_low_0.015.dat"

files = []
for folder in folders_for_infinity_fit:
    print("ON FOLDER:", folder)

    file_list = glob.glob(data_folder + f"{folder[0]}/" + filename)

    for file_path in file_list:
        files.append(file_path)

def analyse_file(file_path):
    print("On file:", file_path)
    data_fit_to_infinity_from_file(file_path, mass=folder[1], dest = file_path)

pool_obj = multiprocessing.Pool()
pool_obj.map(analyse_file, files)
pool_obj.close()
pool_obj.join()

