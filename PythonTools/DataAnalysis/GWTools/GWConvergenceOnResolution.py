
"""
Analyse GW Convergence
(converge of real part, converge of amplitude, convergence of phase, ...)
"""

import sys
sys.path.append("../Base/")

import SmallDataIOReader, IntegrationMethod, DataIntegration
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import make_interp_spline # to interpolate the amplitude
from scipy import integrate # to integrate GW amplitude to align different resolutions
from scipy.optimize import minimize # to align time shift of GW amplitude
from scipy.optimize import fsolve # to find convergence factor
from functools import partial

saveFigs = True
# saveFigs = False
startTime = 300

matplotlib.rcParams.update({'font.size': 18})

# Get Weyl from some file
def get_weyl4_mode(file_path):
    file = SmallDataIOReader.File(file_path)
    columns = file.transposeData(file.getData())
    times = columns[0]
    # keep only first radius (r=50)
    # weyl4_re = columns[1]
    # weyl4_im = columns[2]
    # PICK YOUR RADIUS TO USE
    weyl4_re = columns[11]
    weyl4_im = columns[12]
    weyl4 = [(re, im) for re, im in zip(weyl4_re, weyl4_im)]

    return times, weyl4, weyl4_re, weyl4_im

# Smart conversion to phase such that phase is as continuous as possible
def convert_weyl4_to_amplitude_and_phase(times, weyl4, start_monoticity):
    amplitudes = []
    phases = []
    num_cycles = 0

    total_files = len(weyl4)

    sign = 0

    for i, time, file in zip(range(len(weyl4)), times, weyl4):
        # print("Computation = %d / %d" % (i+1, total_files))

        if isinstance(file, list) and isinstance(file[0], list):
            block0 = file[0]
            weyl4 = block0[1]
        else:
            weyl4 = file

        amplitude = math.sqrt(weyl4[0] ** 2 + weyl4[1] ** 2)
        phase = math.atan2(weyl4[1], weyl4[0])

        # it's pretty tricky to keep the phase monotonous
        # Probably there is a simpler way to do this
        if len(phases) > 0 and time > start_monoticity:
            phase_old = phases[-1]
            phase_minus = phase + (num_cycles-1) * 2 * math.pi
            phase_orig  = phase + (num_cycles) * 2 * math.pi
            phase_plus  = phase + (num_cycles+1) * 2 * math.pi
            diff_minus = abs(phase_minus - phase_old)
            diff_orig  = abs(phase_orig - phase_old)
            diff_plus  = abs(phase_plus - phase_old)

            if diff_orig < diff_minus and diff_orig < diff_plus:
                phase_new = phase_orig
            elif diff_plus < diff_minus:
                num_cycles += 1
                phase_new = phase_plus
            else:
                num_cycles -= 1
                phase_new = phase_minus
        else:
            phase_new = phase

        # print(weyl4, amplitude, phase)
        amplitudes.append(amplitude)
        phases.append(phase_new)

    return amplitudes, phases

def get_amplitude_and_phase_interpolated(times, weyl4, start_monoticity):
    amplitudes, phases = convert_weyl4_to_amplitude_and_phase(times, weyl4, start_monoticity)

    amp_interp = make_interp_spline(times, amplitudes)
    phase_interp = make_interp_spline(times, phases)

    return amplitudes, phases, amp_interp, phase_interp

# Nice function to compute:
# (int of amp1*amp2)/sqrt((int of amp1^2)(int of amp2^2))
# When trying to align functions, we try to minimize this
# This is useful for convergence, this is also useful for certain extrapolations
def amplitude_alignment_factor(amp1, amp2, min_t, max_t, offset):
    if offset < 0:
        min_t_new = min_t - offset
        max_t_new = max_t
    else:
        min_t_new = min_t
        max_t_new = max_t - offset

    assert min_t_new >= min_t
    assert max_t_new <= max_t

    amp1_times_amp2 = integrate.quad(lambda t : amp1(t) * amp2(t + offset), min_t_new, max_t_new)
    amp1_squared = integrate.quad(lambda t : amp1(t) **2, min_t_new, max_t_new)
    amp2_squared = integrate.quad(lambda t : amp2(t + offset) **2, min_t_new, max_t_new)

    # require less than 0.1% error in the integral
    # print(offset, amp1_times_amp2, amp1_times_amp2[1] / amp1_times_amp2[0])
    # assert amp1_times_amp2[1] / amp1_times_amp2[0] < 1e-3
    # assert amp1_squared[1] / amp1_squared[0] < 1e-3
    # assert amp2_squared[1] / amp2_squared[0] < 1e-3

    return amp1_times_amp2[0] / math.sqrt(amp1_squared[0] * amp2_squared[0])

###############
# USING Weyl4ExtractionOut_*.dat

# file_prefix_low = "run_g2_v1_h096/Weyl4ExtractionOut_*.dat"
# # file_prefix_low = "run_g2_v1_h096/Weyl4ExtractionOut_00000?.dat"
# file_prefix_med = "run_g2_v2_h112/Weyl4ExtractionOut_*.dat"

# Not in use anymore
# def get_weyl4_space_integrated(file_prefix):

#     reader = SmallDataIOReader.FileSet(file_prefix, read = False, skip = 40)

#     total_files = reader.numFiles()
#     print("Number of files =", total_files)

#     times = []

#     # keep only first block (r=50)
#     for i, file in enumerate(reader):
#         print("File = %d / %d" % (i+1, total_files))
#         file.read()
#         file.removeBlock(1, file.numBlocks()-1)
#         # print(file.numBlocks())
#         time = file.getHeaderValue("time")
#         # print(time)
#         times.append(time)

#     print("Integrating data in space")
#     weyl4_space_integrated = DataIntegration.integrate_in_space_2d(reader,
#                                             [IntegrationMethod.simpson, IntegrationMethod.trapezium],
#                                             [False, True],
#                                             'r',
#                                             DataIntegration.spherical_area_element)
#     print("Done") 

#     # print(weyl4_space_integrated)
#     # print(type(weyl4_space_integrated))

#     return times, weyl4_space_integrated

# times_low, weyl4_low = get_weyl4_space_integrated(file_prefix_low)
# times_med, weyl4_med = get_weyl4_space_integrated(file_prefix_med)

###############
# USING Weyl_integral_22.dat

Ns = np.array([96, 112, 128])

base_folder = "./run17/"
file_prefix_low = base_folder + "run_g2_v1_h096/weyl4/Weyl_integral_22.dat"
file_prefix_med = base_folder + "run_g2_v2_h112/weyl4/Weyl_integral_22.dat"
file_prefix_high = base_folder + "run_g2_v3_h128/weyl4/Weyl_integral_22.dat"

times_low, weyl4_low, weyl4_re_low, weyl4_im_low = get_weyl4_mode(file_prefix_low)
times_med, weyl4_med, weyl4_re_med, weyl4_im_med = get_weyl4_mode(file_prefix_med)
times_high, weyl4_high, weyl4_re_high, weyl4_im_high = get_weyl4_mode(file_prefix_high)

assert times_low[0] == times_high[0]
assert times_med[0] == times_high[0]
# time_min = max(max(times_low[0] , times_med[0]), times_high[0])
time_min = times_low[0]
time_max = min(min(times_low[-1] , times_med[-1]), times_high[-1])
print("MIN TIME =", time_min)
print("MAX TIME =", time_max)

# This supposedly returns the converge order, but it's super noisy
# def getConvergenceOrder(datas):
#     ratio = abs((datas[0] - datas[1])/(datas[1] - datas[2]))
#     initial_guess = 3.14 # just because (like this I can detect that this is not a coincidence)
#     return fsolve(func, initial_guess, args=(ratio))[0]
# convergenceOrder = [getConvergenceOrder([low, med, high]) for low, med, high in zip(amp_low_offset, amp_med_offset, amp_high_new)]
# print(convergenceOrder)

# Get convergence factor based on 'Ns' (3 grid resolutions)
func = lambda n, c : (np.power(1./Ns[0],n)-np.power(1./Ns[1],n))/\
                        (np.power(1./Ns[1],n)-np.power(1./Ns[2],n)) - c
q2_factor = func(2, 0)
q4_factor = func(4, 0)

###############

def weyl_re_analysis():
    weyl4_re_low_interp = make_interp_spline(times_low, weyl4_re_low)
    weyl4_re_med_interp = make_interp_spline(times_med, weyl4_re_med)
    weyl4_re_high_interp = make_interp_spline(times_high, weyl4_re_high)

    function_to_maximize_low_high = partial(amplitude_alignment_factor, weyl4_re_high_interp, weyl4_re_low_interp, time_min, time_max)
    function_to_maximize_med_high = partial(amplitude_alignment_factor, weyl4_re_high_interp, weyl4_re_med_interp, time_min, time_max)

    minimum_low_high = minimize(lambda offset : -function_to_maximize_low_high(offset), -8)
    minimum_med_high = minimize(lambda offset : -function_to_maximize_med_high(offset), -3)

    best_offset_low_high = minimum_low_high.x
    best_offset_med_high = minimum_med_high.x

    # best_offset_med_high = 0
    # best_offset_low_high = 0

    print("Offset LOW-HIGH =", best_offset_low_high)
    print("Offset MED-HIGH =", best_offset_med_high)

    times_new = np.arange(time_min, time_max, time_max/10000)
    weyl4_re_high_new = weyl4_re_high_interp(times_new)
    weyl4_re_low_offset = weyl4_re_low_interp(times_new + best_offset_low_high)
    weyl4_re_med_offset = weyl4_re_med_interp(times_new + best_offset_med_high)

    weyl4_re_low_med = [abs(low - med) for low, med in zip(weyl4_re_low_offset, weyl4_re_med_offset)]
    weyl4_re_med_high = [abs(med - high) for med, high in zip(weyl4_re_med_offset, weyl4_re_high_new)]
    q2_weyl4_re_high = [ high * q2_factor for high in weyl4_re_med_high]
    q4_weyl4_re_high = [ high * q4_factor for high in weyl4_re_med_high]

    ##################################################################
    # PLOT Weyl4 Real part
    fig1, ax1 = plt.subplots(figsize=(16, 9))

    ax1.plot(times_new, weyl4_re_high_new, label = "High")
    ax1.plot(times_new, weyl4_re_low_offset, "--", label = "Low Shifted")
    ax1.plot(times_new, weyl4_re_med_offset, "--", label = "Medium Shifted")
    # fig1.suptitle('GW Weyl4 Real Part', position=(0.5,0.93))

    # ax1.set_yscale("log")

    ax1.set_xlabel(r'$t~(M)$')
    # ax1.ylabel(r"Weyl$_4$ mode-22 Real Part")

    ax1.legend(loc="upper left")

    if saveFigs:
        plt.draw()
        out_name = 'GW_weyl4_re.png'
        plt.savefig(out_name, bbox_inches = 'tight')

    plt.show()

    ##################################################################
    # PLOT Weyl4 Real Part CONVERGENCE
    fig2, ax2 = plt.subplots(figsize=(16, 9))

    ax2.plot(times_new, weyl4_re_low_med, label = "LOW-MED")
    ax2.plot(times_new, weyl4_re_med_high, label = "MED-HIGH")
    ax2.plot(times_new, q2_weyl4_re_high, "--", label = "Q2(MED-HIGH)")
    ax2.plot(times_new, q4_weyl4_re_high, "--", label = "Q4(MED-HIGH)")
    # fig2.suptitle('GW Weyl4 Real Part Convergence', position=(0.5,0.93))

    ax2.set_yscale("log")

    ax2.set_xlabel(r'$t~(M)$')
    # ax2.ylabel(r"Weyl$_4$ mode-22 Real Part")

    ax2.legend(loc="upper left")

    ax2ins = ax2.inset_axes([0.4, 0.08, 0.5, 0.4])
    ax2ins.plot(times_new, weyl4_re_low_med, label = "LOW-MED")
    ax2ins.plot(times_new, weyl4_re_med_high, label = "MED-HIGH")
    ax2ins.plot(times_new, q2_weyl4_re_high, "--", label = "Q2(MED-HIGH)")
    ax2ins.plot(times_new, q4_weyl4_re_high, "--", label = "Q4(MED-HIGH)")
    ax2ins.set_xlim(1225, 1325)
    ax2ins.set_ylim(1e-5, 5e-3)
    ax2ins.set_yscale("log")
    # ax2ins.set_xticklabels('')
    # ax2ins.set_yticklabels('')
    ax2.indicate_inset_zoom(ax2ins, edgecolor="black")

    if saveFigs:
        plt.draw()
        out_name = 'GW_weyl4_re_convergece.png'
        plt.savefig(out_name, bbox_inches = 'tight')

    plt.show()


def amplitude_and_phase_analysis():
    amp_low, phase_low, amp_low_interp, phase_low_interp = get_amplitude_and_phase_interpolated(times_low, weyl4_low, startTime)
    amp_med, phase_med, amp_med_interp, phase_med_interp = get_amplitude_and_phase_interpolated(times_med, weyl4_med, startTime)
    amp_high, phase_high, amp_high_interp, phase_high_interp = get_amplitude_and_phase_interpolated(times_high, weyl4_high, startTime)

    function_to_maximize_low_high = partial(amplitude_alignment_factor, amp_high_interp, amp_low_interp, time_min, time_max)
    function_to_maximize_med_high = partial(amplitude_alignment_factor, amp_high_interp, amp_med_interp, time_min, time_max)

    # print(function_to_maximize_low_high(0))
    # print(function_to_maximize_low_high(0.1))
    # print(function_to_maximize_low_high(-0.1))

    minimum_low_high = minimize(lambda offset : -function_to_maximize_low_high(offset), -8)
    minimum_med_high = minimize(lambda offset : -function_to_maximize_med_high(offset), -3)

    best_offset_low_high = minimum_low_high.x
    best_offset_med_high = minimum_med_high.x

    print("Offset LOW-HIGH =", best_offset_low_high)
    print("Offset MED-HIGH =", best_offset_med_high)

    times_new = np.arange(time_min, time_max, time_max/10000)
    # amp_low_new = amp_low_interp(times_new)
    # amp_med_new = amp_med_interp(times_new)
    amp_high_new = amp_high_interp(times_new)
    amp_low_offset = amp_low_interp(times_new + best_offset_low_high)
    amp_med_offset = amp_med_interp(times_new + best_offset_med_high)

    # phase_low_new = phase_low_interp(times_new)
    # phase_med_new = phase_med_interp(times_new)
    phase_high_new = phase_high_interp(times_new)
    phase_low_offset = phase_low_interp(times_new + best_offset_low_high)
    phase_med_offset = phase_med_interp(times_new + best_offset_med_high)
    # phase_low_offset = phase_low_interp(times_new)
    # phase_med_offset = phase_med_interp(times_new)

    amp_low_med = [abs(low - med) for low, med in zip(amp_low_offset, amp_med_offset)]
    amp_med_high = [abs(med - high) for med, high in zip(amp_med_offset, amp_high_new)]
    q2_amp_high = [ high * q2_factor for high in amp_med_high]
    q4_amp_high = [ high * q4_factor for high in amp_med_high]

    phase_low_med = [abs(low - med) for low, med in zip(phase_low_offset, phase_med_offset)]
    phase_med_high = [abs(med - high) for med, high in zip(phase_med_offset, phase_high_new)]
    q2_phase_high = [ high * q2_factor for high in phase_med_high]
    q4_phase_high = [ high * q4_factor for high in phase_med_high]

    ##################################################################
    # PLOT AMPLITUDE
    fig1, ax1 = plt.subplots(figsize=(16*2, 9*2))

    # ax1.plot(times_new, amp_low_new, label = "Low")
    # ax1.plot(times_new, amp_med_new, label = "Medium")
    ax1.plot(times_new, amp_high_new, label = "High")
    ax1.plot(times_new, amp_low_offset, "--", label = "Low Shifted")
    ax1.plot(times_new, amp_med_offset, "--", label = "Medium Shifted")
    # ax1.plot(times_low, amp_low, label = "Low")
    # ax1.plot(times_med, amp_med, label = "High")
    fig1.suptitle('GW Amplitude', position=(0.5,0.93))

    ax1.set_yscale("log")

    ax1.set_xlabel(r'$t~(M)$')
    # ax1.ylabel(r"Weyl$_4$ mode-22 Amplitude")

    ax1.legend()

    if saveFigs:
        plt.draw()
        out_name = 'GW_amplitude.png'
        plt.savefig(out_name, bbox_inches = 'tight')

    plt.show()

    ##################################################################
    # PLOT AMPLITUDE CONVERGENCE
    fig2, ax2 = plt.subplots(figsize=(16*2, 9*2))

    ax2.plot(times_new, amp_low_med, label = "LOW-MED")
    ax2.plot(times_new, amp_med_high, label = "MED-HIGH")
    ax2.plot(times_new, q2_amp_high, "--", label = "Q2(MED-HIGH)")
    ax2.plot(times_new, q4_amp_high, "--", label = "Q4(MED-HIGH)")
    fig2.suptitle('GW Amplitude Convergence', position=(0.5,0.93))

    ax2.set_yscale("log")

    ax2.set_xlabel(r'$t~(M)$')
    # ax2.ylabel(r"Weyl$_4$ mode-22 Amplitude")

    ax2.legend()

    if saveFigs:
        plt.draw()
        out_name = 'GW_amplitude_convergece.png'
        plt.savefig(out_name, bbox_inches = 'tight')

    plt.show()

    ##################################################################
    # PLOT PHASE
    fig3, ax3 = plt.subplots(figsize=(16*2, 9*2))

    # ax3.plot(times_low, phase_low, label = "LOW")
    # ax3.plot(times_med, phase_med, label = "MED")
    ax3.plot(times_new, phase_low_offset, label = "LOW")
    ax3.plot(times_new, phase_med_offset, label = "MED")
    ax3.plot(times_new, phase_high_new, label = "HIGH")
    fig3.suptitle('GW Phase', position=(0.5,0.93))

    # ax3.yscale("log")

    ax3.set_xlabel(r'$t~(M)$')
    # ax3.ylabel(r"Weyl$_4$ mode-22 Phase")

    ax3.legend()

    if saveFigs:
        plt.draw()
        out_name = 'GW_phase.png'
        plt.savefig(out_name, bbox_inches = 'tight')

    plt.show()

    ##################################################################
    # PLOT PHASE CONVERGENCE
    fig4, ax4 = plt.subplots(figsize=(16*2, 9*2))

    ax4.plot(times_new, phase_low_med, label = "LOW-MED")
    ax4.plot(times_new, phase_med_high, label = "MED-HIGH")
    ax4.plot(times_new, q2_phase_high, "--", label = "Q2(MED-HIGH)")
    ax4.plot(times_new, q4_phase_high, "--", label = "Q4(MED-HIGH)")
    fig4.suptitle('GW Phase Convergence', position=(0.5,0.93))

    ax4.set_yscale("log")

    ax4.set_xlabel(r'$t~(M)$')
    # ax4.ylabel(r"Weyl$_4$ mode-22 Phase")

    ax4.legend()

    if saveFigs:
        plt.draw()
        out_name = 'GW_phase_convergence.png'
        plt.savefig(out_name, bbox_inches = 'tight')

    plt.show()

weyl_re_analysis()
# amplitude_and_phase_analysis()
