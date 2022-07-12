
from scipy import integrate # to integrate GW amplitude to align different resolutions
import math
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from scipy.optimize import minimize # to align time shift of GW amplitude
from scipy.interpolate import make_interp_spline # to interpolate the amplitude

import sys
sys.path.append("./Base/")
import SmallDataIOReader

################################################################################
# Setting the font structure to latex-type
# Just call 'use_tex' before printing to pdf with matplotlib
from pylab import *
resize = 2
def use_tex(fontsize=44, legend_fontsize = 32):
    rc = mpl.rcParams # Font structure is called rc now
    rc['text.usetex'] = True # Tex fonts
    rc['font.family'] = 'sans-serif'
    rc['font.serif'].insert(0,'cm') # Default font is computer modern for latex
    rc['font.size'] = fontsize*resize # Axes font size
    rc['xtick.labelsize'] = 'small'
    rc['ytick.labelsize'] = 'small'
    rc['legend.fontsize'] = legend_fontsize*resize # Legend font size
################################################################################

def tortoise_radius(radius, mass): return radius + 2*mass*math.log(radius/(2*mass)-1)

# Get the headers radius from .dat GRChombo files
# (i.e. the "r = 50  50  60  60  ..." in our .dat files)
def get_radii(file, verbose=True, has_infinity=False):
    f = open(file,'r')
    headers = f.readline()
    radii = f.readline()
    f.close()

    headers = headers.split()[2:] # remove comment and 'time'
    headers = list(dict.fromkeys(headers))
    n_vars = len(headers)
    if verbose:
        print("Headers: ", headers)

    radii = radii.split()[3:] # remove comment and 'r' and '='
    radii = list(dict.fromkeys(radii))
    if has_infinity:
        radii = radii[:-1] # remove Infinity
    radii = [float(radius) for radius in radii]
    n_radii = len(radii)
    if verbose:
        print("Radii: ", radii)

    return radii

# Richardson of 2nd or 3rd order
# 2nd order if length 2, 3rd order if more points/lists are passed
def richardson_extrapolation(rs, data):
    if len(data) > 2:
        c2 = rs[-2] / rs[-1] * (1. - rs[-3] / rs[-1]) / (1. - rs[-3] / rs[-2]);
        c3 = -rs[-3] / rs[-1] * (1. - rs[-2] / rs[-1]) / (rs[-2] / rs[-3] - 1.);
        comp_inf = (data[-1] - c2 * data[-2] - c3 * data[-3]) / (1. - c2 - c3);
    else:
        comp_inf = (data[-1] - rs[-2] / rs[-1] * data[-2]) / (1. - rs[-2] / rs[-1]);
    return comp_inf

"""
Generic fitting function
E.g.
xs = [1,2,3,4,5]
ys = [1,4,9,16,25]
def func(params, xs):
    a,b,c = params
    return a + b * xs + c * (xs**2)
n_params = 3
guesses_abc = [1,1,1]
params, errors = FIT(xs,ys, func, n_params, beta0 = guesses_abc)
"""
def FIT(x, y, func, nparams, err_x=None, err_y=None, beta0=0, verbose=True):
    # Model object
    from scipy import odr
    model = odr.Model(func)

    if err_x is None and err_y is None:
        data = odr.RealData(x, y)
    elif err_x is None:
        data = odr.RealData(x, y, sy=err_y)
    elif err_y is None:
        data = odr.RealData(x, y, sx=err_x)
    else:
        data = odr.RealData(x, y, sx=err_x, sy=err_y)

    # Set up ODR with the model and data.
    if beta0 == 0:
        beta0 = [1.]*nparams
    odr = odr.ODR(data, model, beta0=beta0)
    out = odr.run()

    #print fit parameters and 1-sigma estimates
    popt = out.beta
    perr = out.sd_beta
    if verbose:
        print('Fit Parameters 1-sigma error')
        for i in range(len(popt)):
            print(str(popt[i])+' +- '+str(perr[i]))
    return popt, perr


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

# Nice function to align waves. Used for convergence tests and certain extrapolations
def align_waves(data_base, datas_other, time_base, times_other, initial_guesses, pair_base=None, pairs_other=None):
    data_base_interp = make_interp_spline(time_base, data_base)
    datas_other_interp = [make_interp_spline(time, data) for time, data in zip(times_other, datas_other)]

    if pair_base is not None:
        pair_base_interp = make_interp_spline(time_base, pair_base)
        pairs_other_interp = [make_interp_spline(time, data) for time, data in zip(times_other, pairs_other)]

    time_min = max(np.max([np.min(time) for time in times_other]), np.min(time_base))
    time_max = min(np.min([np.max(time) for time in times_other]), np.max(time_base))

    # print("time_min =", time_min)
    # print("time_max =", time_max)

    functions_to_maximize = [partial(amplitude_alignment_factor, data_base_interp, data_other_interp, time_min, time_max)
                                        for data_other_interp in datas_other_interp]
    minima = [minimize(lambda offset : -function_to_maximize(offset), -8) for function_to_maximize in functions_to_maximize] 

    best_offsets = [minimum.x for minimum in minima]

    print("Offsets =", best_offsets)

    time_new = np.linspace(time_min, time_max, 10000) 
    data_base_new = data_base_interp(time_new)
    datas_other_new = [data_other_interp(time_new + best_offset) for data_other_interp, best_offset in zip(datas_other_interp, best_offsets)]

    if pairs_other is not None:
        pair_base_new = pair_base_interp(time_new)
        pairs_other_new = [pair_other_interp(time_new + best_offset) for pair_other_interp, best_offset in zip(pairs_other_interp, best_offsets)]
        return time_new, np.append(datas_other_new, [data_base_new], axis=0), np.append(pairs_other_new, [pair_base_new], axis=0), best_offsets

    else:
        return time_new, np.append(datas_other_new, [data_base_new], axis=0), best_offsets


"""
Some helper functions for discrete fourier transform
"""
def zero_pad(data):
    """
    This function takes in a vector and zero pads it so it is a power of two.
    """
    N = len(data)
    # use A LOT of padding, makes FFT smoother (smaller frequency 'df')
    # (recall that 'df' of FFT is proportional to total time of the data
    # so we pad with zeros such that the FFT has more resolution)
    pow_2 = np.ceil(np.log2(N))+3
    padded = np.pad(data,(0,int((2**pow_2)-N)),'constant')
    return padded

def FFT(waveform):
    """
    Here we taper the signal, pad and then compute the FFT.
    """
    N = len(waveform)
    # taper with tukey window to ensure that extremes of the wave go to 0
    # (https://en.wikipedia.org/wiki/Window_function#Tukey_window)
    taper = scipy.signal.tukey(N, 0.04) # ~50M for a GW of ~1400M
    # now add zeros to increase FFT resolution (see comments inside zero_pad function)
    waveform_w_pad = zero_pad(waveform * taper)
    roll = waveform_w_pad.size // 2
    # we need to 'roll' the wave because 'fft' default returns the frequencies rotated
    # (from 0 to 1 instead of -1/2 to 1/2 for negative and positive frequencies)
    frequencies = np.roll( np.fft.fftfreq(waveform_w_pad.size) , roll )
    fft = np.roll( np.fft.fft(waveform_w_pad) , roll )
    return frequencies, fft

def iFFT(waveform_w):
    roll = waveform_w.size // 2
    return np.fft.ifft(np.roll(waveform_w, -roll))


def plot_FFT_from_file(filepaths, column = 1, use_dts=False, max_frequency=None):
    waves = []
    dts = []
    for filepath in filepaths:
        mode = SmallDataIOReader.File(filepath).getBlock(0)
        mode_data = np.transpose(mode.getData())
        waves.append(mode_data[column] + 1j * mode_data[column+1])
        dts.append((mode_data[0][1] - mode_data[0][0]) if use_dts else 1)
    plot_FFT(waves, dts, max_frequency=max_frequency)


def plot_FFT(waves, dts=None, max_frequency = None):
    if dts is None:
        dts = [1] * len(waves)
    for i, dt, wave in zip(range(len(waves)), dts, waves):
        freqs, mode_FFT = FFT(wave)
        freqs /= dt
        # plt.plot(freqs, np.real(mode_FFT) / len(wave), label=f"Re Wave {i}")
        plt.plot(freqs, np.abs(mode_FFT) / len(freqs), label=f"Re Wave {i}")
        # plt.plot(freqs, np.imag(mode_FFT) / len(wave), label=f"Im Wave {i}")
    
    if len(waves) > 1:
        plt.legend()
    if max_frequency is not None:
        plt.xlim(-max_frequency, max_frequency)
    plt.show()


def convert_complex_to_amplitude_and_phase(times, values, start_monoticity):
    num_cycles = 0
    # sign = 0

    amplitudes = np.abs(values)
    phases = np.angle(values)

    for i, time, value in zip(range(len(values)), times, values):
        if i > 0 and time >= start_monoticity:
            phase_old = phases[i-1]
            phase = phases[i]

            phase_minus = phase + (num_cycles-1) * 2 * math.pi
            phase_orig  = phase + (num_cycles) * 2 * math.pi
            phase_plus  = phase + (num_cycles+1) * 2 * math.pi
            diff_minus = abs(phase_minus - phase_old)
            diff_orig  = abs(phase_orig - phase_old)
            diff_plus  = abs(phase_plus - phase_old)

            if diff_orig < diff_minus and diff_orig < diff_plus:
                phase_new = phase_orig
            elif diff_plus < diff_minus:
                # if sign != -1: # always cycle in the same direction
                # sign = 1
                num_cycles += 1
                phase_new = phase_plus
            else:
                # if sign != 1: # always cycle in the same direction
                # sign = -1
                num_cycles -= 1
                phase_new = phase_minus

            # print(time, phase_minus, phase_orig, phase_plus, phase_new, num_cycles)
            
            phases[i] = phase_new

    return amplitudes, phases

