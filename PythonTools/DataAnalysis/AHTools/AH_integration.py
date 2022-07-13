
"""
Example file integrating quantities over the AH coords files
"""

import sys
sys.path.append("../Base/")

import SmallDataIOReader, IntegrationMethod, DataIntegration, SphericalHarmonics

import numpy as np
import math
import matplotlib.pyplot as plt
import os

labels = ["g2", "GR"]
location = ["../run_g2_v3_h128/ah/", "../run_GR_v3_h128/ah/"]

collects = ["phi", "WeakField", "rho"]
names = ["Phi", "Weak Field Condition", "Energy density"]


#####################################################################

def amplitude(x):
    return math.sqrt(x.real ** 2 + x.imag ** 2)

def get_integrated_modes(location, ah_index=1, modes = [(2,2)], jump=1):
    file_prefix = location + "coords_AH" + str(ah_index) + "_*.dat"
    stats_file = location + "stats_AH" + str(ah_index) + ".dat"

    stats = SmallDataIOReader.File(stats_file)

    data = stats.getData()
    data_no_nan = data[~np.isnan(data).any(axis=1)]

    times = data_no_nan[:,0]
    files = data_no_nan[:,1]

    # times = times[:20]
    # files = files[:20]

    total_collects = len(collects)

    total_files = len(times)
    integrated_lm = []
    for c in range(total_collects):
        integrated_lm.append([])
        for mode in modes:
            integrated_lm[c].append([])

    valid_times = []

    num_modes = len(modes)

    for i, time, file_i in zip(range(total_files), times, files):
        if i%jump == 0:
            try:
                print("File = %d / %d at time %f" % (i+1, total_files, time))
                filename = location + "coords_AH" + str(ah_index) + "_%06d.dat" % file_i
                # print(filename)
                file = SmallDataIOReader.File(filename)

                for collect in collects:
                    for mode in modes:
                        file.addColumn(lambda row, map: row[map[collect][1]] * SphericalHarmonics.spin_Y_lm(theta=row[0], phi=row[1], es=0, em=mode[1], el=mode[0]))
                file.addColumn(lambda row, map: 1)
                file.removeColumn(file.numLabels(), (file.numColumns()-1) - total_collects * num_modes - 1)

                integrated_lm_i = DataIntegration.integrate_in_space_2d(file,
                                                                        # [IntegrationMethod.simpson, IntegrationMethod.trapezium],
                                                                        [IntegrationMethod.trapezium, IntegrationMethod.trapezium],
                                                                        [False, True],
                                                                        'r',
                                                                        DataIntegration.spherical_area_element)

                normalization = integrated_lm_i[0][-1]
                assert normalization.imag == 0
                normalization = normalization.real

                for c in range(total_collects):
                    for m in range(num_modes):
                        value = amplitude(integrated_lm_i[0][num_modes * c + m])
                        integrated_lm[c][m].append(value / normalization)

                valid_times.append(time)

            except:
                print("File = %d / %d at time %f FAILED" % (i+1, total_files, time))

    # print(integrated_lm)
    # print(len(integrated_lm))

    return (valid_times, integrated_lm)

modes = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2), (4,0), (4,2), (4,4), (6,0), (6,2), (6,4), (6,6)]
jump = 1 # jump every X files

integrated_phi_lm_ah1 = []
integrated_phi_lm_merger = []
for loc in location:
    integrated_phi_lm_ah1.append( get_integrated_modes(loc, modes=modes, jump=jump) )
    integrated_phi_lm_merger.append( get_integrated_modes(loc, ah_index=3, modes=modes, jump=jump) )

##################################################################
# PLOT
def plot(collect_index, mode_index, saveFigs=False):
    fig1, ax1 = plt.subplots(figsize=(16, 9))

    for l in range(len(location)):
        # for m, mode in enumerate(modes):
        ax1.plot(integrated_phi_lm_ah1[l][0], integrated_phi_lm_ah1[l][1][collect_index][mode_index],       
            color="C" + str(l),
            label = f"l={modes[mode_index][0]}, m={modes[mode_index][1]} ({labels[l]})")
        ax1.plot(integrated_phi_lm_merger[l][0], integrated_phi_lm_merger[l][1][collect_index][mode_index],
            color="C" + str(l))
    fig1.suptitle(f'{names[collect_index]} comparison on AH - mode l={modes[mode_index][0]}, m={modes[mode_index][1]}', position=(0.5,0.93))

    ax1.set_yscale("log")

    ax1.set_xlabel(r'$t~(M)$')
    # ax1.ylabel(r"Weyl$_4$ mode-22 Real Part")

    ax1.legend()

    # normal plot
    folder = f"./{collects[collect_index]}/"
    if saveFigs:
        plt.draw()
        if not os.path.isdir(folder):
            os.mkdir(folder)
        out_name = folder + f"{collects[collect_index]}_modes_l{modes[mode_index][0]}_m{modes[mode_index][1]}.png"
        plt.savefig(out_name, bbox_inches = 'tight')

    # zoomed in plot
    ax1.set_xlim(850, 1350)

    if saveFigs:
        plt.draw()
        out_name = folder + f"{collects[collect_index]}_modes_l{modes[mode_index][0]}_m{modes[mode_index][1]}_zoom.png"
        plt.savefig(out_name, bbox_inches = 'tight')

    # plt.show()
    plt.close()

saveFigs = True
# saveFigs = False

for c in range(len(collects)):
    for m in range(len(modes)):
        plot(collect_index = c, mode_index = m, saveFigs = saveFigs)
