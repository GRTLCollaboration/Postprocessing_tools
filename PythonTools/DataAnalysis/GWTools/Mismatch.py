
import sys, os
sys.path.append("../")
from utils import *
from Weyl4ToStrain import *

from scipy.interpolate import make_interp_spline
from scipy import integrate # to integrate GW amplitude to align different resolutions
from functools import partial
from scipy.optimize import minimize # to align time shift for mismatch
import multiprocessing
import numpy as np

class LIGO_PSD:
    def __init__(self):
        # Updated Advanced LIGO sensitivity design curve
        # from: https://dcc.ligo.org/LIGO-T1800044/public
        ligo_filename = "aLIGODesign.txt"
        freqs, psds = np.loadtxt(ligo_filename, unpack=True)

        self.func = make_interp_spline(freqs, psds, k=3)
        self.minimum = freqs[0]
        self.maximum = freqs[-1]

    def get(self, freq):
        assert freq >= self.minimum
        assert freq <= self.maximum
        return self.func(freq)

class Constant_PSD:
    def __init__(self):
        self.minimum = 0
        self.maximum = 10e10

    def get(self, freq):
        return 1.

G = 6.67e-11 # kg^-1 m^3 s^-2
c = 3e8 # m s^-1
mass_sun = 2e30 # kg
mass_sum_length = G * mass_sun / c**2 # m
mass_sun_time = mass_sum_length / c # s

# from 1006.1632 eqs. 28-30
# def fake_gw(t):
#     amp = 0.02 / 2 * (1 + np.tanh((t+480)/10)) * (1 + 5 * np.exp(t/80) * (1+np.tanh(-(t-0)/16)))
#     phase = 0.2*(t-0) + (1-0.2)/2*(1+80*np.log(np.cosh((t-0)/80)))
#     return amp * np.sin(-phase)

def complex_integral(func, x_min, x_max):
    integral_real, error_real = integrate.quad(lambda x: func(x).real, x_min, x_max)
    integral_imag, error_imag = integrate.quad(lambda x: func(x).imag, x_min, x_max)
    return (integral_real + 1j*integral_imag, error_real + 1j*error_imag)

def inner_product(time_min, time_max, dt, wave_1_re, wave_1_im, wave_2_re, wave_2_im,
    solar_masses = 100, cut_off_high_freq = 5000, verbose=True,
    PSD = LIGO_PSD()):

    # use whatever you want, atm uncommented below is the real part, i.e. we are
    # computing the mismatch of the + polarization, h+

    # wave 1
    freqs_1, wave_1_re_FFT = FFT(wave_1_re)
    # freqs_1_im, wave_1_im_FFT = FFT(wave_1_im)
    # assert len(freqs_1) == len(freqs_1_im)

    freqs_1 /= dt
    # wave_1_FFT_interp = make_interp_spline(freqs_1, wave_1_re_FFT + 1j * wave_1_im_FFT, k=3)
    wave_1_FFT_interp = make_interp_spline(freqs_1, wave_1_re_FFT, k=3)

    # wave 2
    freqs_2, wave_2_re_FFT = FFT(wave_2_re)
    # freqs_2_im, wave_2_im_FFT = FFT(wave_2_im)
    # assert len(freqs_2) == len(freqs_2_im)

    freqs_2 /= dt
    # wave_2_FFT_interp = make_interp_spline(freqs_2, wave_2_re_FFT + 1j * wave_2_im_FFT, k=3)
    wave_2_FFT_interp = make_interp_spline(freqs_2, wave_2_re_FFT, k=3)

    # inner product
    freq_min = max(freqs_1[1]-freqs_1[0], freqs_2[1]-freqs_2[0])
    freq_max = min(min(freqs_1[-1], freqs_2[-1]), cut_off_high_freq)

    time_convertion = mass_sun_time * solar_masses

    # if verbose:
    #     print("Min/Max frequency:", freq_min / time_convertion, "/", freq_max / time_convertion)
    #     print("PSD min/max:", PSD.minimum, "/", PSD.maximum)

    freq_min = max(freq_min, PSD.minimum * time_convertion)
    freq_max = min(freq_max, PSD.maximum * time_convertion)

    if verbose:
        print("Using min/max:", freq_min / time_convertion, "/", freq_max / time_convertion)

    # plt.plot(freqs_1, wave_1_re_FFT / len(freqs_1))
    # plt.plot(freqs_2, wave_2_re_FFT / len(freqs_2))
    # plt.show()
    # exit()

    integral, error = complex_integral(lambda f : 4 * wave_1_FFT_interp(f).conj() * wave_2_FFT_interp(f) / PSD.get(f / time_convertion),
                                    freq_min, freq_max)

    # computing the absolute value instead of the real part is equivalent to maximizing the real part with a phase added to wave_2
    return np.abs(integral), np.abs(error), time_min, time_max, freq_min, freq_max


def wave_overlap(times_1, wave_1_re, wave_1_im, times_2, wave_2_re, wave_2_im,
    solar_masses = 100, cut_off_high_freq = 5000, verbose=True, PSD = LIGO_PSD(), label="", time_shift = None):
    if time_shift is None:
        time_shift = 0
    elif verbose:
        print("##############################")
        print("Using time_shift =", time_shift, ", solar_masses = ", solar_masses)
    # minimize sends an array, not a single float
    time_shift = time_shift[0] if hasattr(time_shift, '__len__') else time_shift

    # use same time for all
    time_min = max(max(times_1[0], times_2[0] - time_shift), times_2[0])
    time_max = min(min(times_1[-1], times_2[-1] - time_shift), times_2[-1])

    # if verbose:
    #     print(f"Time interval = [{time_min:.4f},{time_max:.4f}]")

    times = np.linspace(time_min, time_max, 10000)
    dt = times[1] - times[0]

    wave_1_re_interp = make_interp_spline(times_1, wave_1_re, k=3)
    wave_1_re_full = wave_1_re_interp(times)
    wave_1_im_interp = make_interp_spline(times_1, wave_1_im, k=3)
    wave_1_im_full = wave_1_im_interp(times)

    wave_2_re_interp = make_interp_spline(times_2, wave_2_re, k=3)
    wave_2_re_full = wave_2_re_interp(times + time_shift)
    wave_2_im_interp = make_interp_spline(times_2, wave_2_im, k=3)
    wave_2_im_full = wave_2_im_interp(times + time_shift)

    inner_product_1_2, error_1_2 = inner_product(time_min, time_max, dt, wave_1_re_full, wave_1_im_full, wave_2_re_full, wave_2_im_full,
                                solar_masses=solar_masses, cut_off_high_freq=cut_off_high_freq, verbose=verbose, PSD=PSD)[:2]
    inner_product_1_1, error_1_1 = inner_product(time_min, time_max, dt, wave_1_re_full, wave_1_im_full, wave_1_re_full, wave_1_im_full,
                                solar_masses=solar_masses, cut_off_high_freq=cut_off_high_freq, verbose=False, PSD=PSD)[:2]
    inner_product_2_2, error_2_2 = inner_product(time_min, time_max, dt, wave_2_re_full, wave_2_im_full, wave_2_re_full, wave_2_im_full,
                                solar_masses=solar_masses, cut_off_high_freq=cut_off_high_freq, verbose=False, PSD=PSD)[:2]
    if verbose:
        print(f"Errors = {error_1_2 / inner_product_1_2 * 100:.2f} % / {error_1_1 / inner_product_1_1 * 100:.2f} % / {error_2_2 / inner_product_2_2 * 100:.2f} %")

    normalized = inner_product_1_2 / sqrt(inner_product_1_1 * inner_product_2_2)
    if normalized > 1:
        # This does happen sometimes. You will see it if it happens.
        print("WARNINGGGGGGGGGGGGGGGGGGGGG")
        print(f"Got {normalized:.4f} from ({inner_product_1_2}+-{error_1_2}), ({inner_product_1_1}+-{error_1_1}), ({inner_product_2_2}+-{error_2_2})")
        print("WARNINGGGGGGGGGGGGGGGGGGGGG")
        with open("mismatch_warnings.dat", 'a') as file:
            file.write(f"Warning for mass={solar_masses}, cutoff={cut_off_high_freq}, label={label}:\n")
            file.write(f"Got {normalized:.4f} from ({inner_product_1_2}+-{error_1_2}), ({inner_product_1_1}+-{error_1_1}), ({inner_product_2_2}+-{error_2_2})\n")
        # assert normalized <= 1
    return normalized, inner_product_1_2, inner_product_1_1, inner_product_2_2

def maximize_wave_overlap(obj):
    times_1, wave_1_re, wave_1_im, times_2, wave_2_re, wave_2_im, solar_masses, cut_off_high_freq, verbose, PSD, label, guess = obj
    print(f"Guess = {guess} (mass = {solar_masses:.2f}, label = {label})")
    function_to_maximize = partial(wave_overlap, times_1, wave_1_re, wave_1_im, times_2, wave_2_re, wave_2_im, solar_masses, cut_off_high_freq, verbose, PSD, label)
    return minimize(lambda time_shift : -function_to_maximize(time_shift)[0], guess)

def compute_mismatch(times_1, wave_1_re, wave_1_im, times_2, wave_2_re, wave_2_im, solar_masses = 100, cut_off_high_freq = 5000, verbose=True, PSD=LIGO_PSD(), parallel=True, label=""):
    obj = [times_1, wave_1_re, wave_1_im, times_2, wave_2_re, wave_2_im, solar_masses, cut_off_high_freq, False, PSD, label]

    # to compute the mismatch, we minimize over time shifts
    # this minimization is often sensitive to the initial guess, so we provide several guesses
    # and then hope one of them convergences to the absolute minimum
    # typically all of them converge to the same

    # you may have to play with the initial guesses for the time shifts if you see that
    # the solver is not converging nicely
    # guesses = [-15, -8, -1, 1, 8, 15]

    # we use the difference in the peaks as initial guess below:    
    wave_1_amp = np.abs(wave_1_re + 1j * wave_1_im)
    wave_2_amp = np.abs(wave_2_re + 1j * wave_2_im)
    time_diff = times_2[np.argmax(wave_2_amp)] - times_1[np.argmax(wave_1_amp)]
    guesses = [time_diff, time_diff * 1.25, time_diff * 0.75, 0]
    print("TIME DIFF =", time_diff)

    if parallel:
        pool_obj = multiprocessing.Pool()
        minima = pool_obj.map(maximize_wave_overlap, [obj + [guess] for guess in guesses])
        pool_obj.close()
        pool_obj.join()
    else:
        minima = [maximize_wave_overlap(obj + [guess]) for guess in guesses]

    print(f"Converged to: [{[(1+minimum.fun)*100 for minimum in minima]}]\n\tat: [{[minimum.x[0] for minimum in minima]}]")
    minimum = min(minima, key=lambda m: m.fun)
    best_offset = minimum.x[0]
    print("Offsets =", best_offset)

    # recompute the optimal one
    normalized, inner_product_1_2, inner_product_1_1, inner_product_2_2 = wave_overlap(times_1, wave_1_re, wave_1_im, times_2, wave_2_re, wave_2_im,
    solar_masses = solar_masses, cut_off_high_freq = cut_off_high_freq, verbose=verbose, PSD = PSD, label=label, time_shift = best_offset)

    # the overlap is -minimum.fun, so mismatch=1-overlap=1+minimum.fun
    return np.array([1 + minimum.fun, inner_product_1_2, inner_product_1_1, inner_product_2_2, best_offset, time_diff])

def compute_mismatch_multiprocessing_wrapper(obj):
    strain_GR, strain_other, cut_off_high_freq, PSD, label, mass = obj
    return compute_mismatch(strain_GR[0], strain_GR[3], strain_GR[4], strain_other[0], strain_other[3], strain_other[4],
            cut_off_high_freq=cut_off_high_freq,
            solar_masses = mass,
            parallel=False,
            PSD=PSD,
            label=label)

def compute_all_mismatch_data(file_GR, file_other, cut_off_high_freq=0.2):
    print("ON FILE:", file_other)

    strain_GR = np.transpose(np.loadtxt(file_GR))
    strain_other = np.transpose(np.loadtxt(file_other))
    label = file_other

    mismatch_constant = compute_mismatch(strain_GR[0], strain_GR[3], strain_GR[4], strain_other[0], strain_other[3], strain_other[4],
                cut_off_high_freq=cut_off_high_freq,
                solar_masses = 100, # doesn't matter
                parallel=True,
                PSD=Constant_PSD(),
                label=label)

    # masses from 1 to ~1000 (2*10=1024)
    # use 10.01 for the np.arrange to include '10'
    solar_masses = np.power(2, np.arange(0, 10.01, 0.5))
    # solar_masses = np.array([100])

    obj = [strain_GR, strain_other, cut_off_high_freq, LIGO_PSD(), label]

    pool_obj = multiprocessing.Pool()
    mismatches = pool_obj.map(compute_mismatch_multiprocessing_wrapper, [obj + [mass] for mass in solar_masses])
    pool_obj.close()
    pool_obj.join()

    solar_masses = np.append(solar_masses, 0)
    mismatches = np.append(mismatches, [mismatch_constant], axis=0)
    mismatches = np.array([np.append([mass], mismatch) for mass, mismatch in zip(solar_masses, mismatches)])

    print("Mismatch constant =", mismatch_constant[0]*100, "%")
    print("Mismatch LIGO =", mismatches[:-1,1]*100, "%")
    print("(masses =", mismatches[:-1,0],")")

    # filename = os.path.dirname(file_other) + "/mismatches_vs_g3_-0.03.dat"
    filename = os.path.dirname(file_other) + "/mismatches_vs_GR.dat"
    # filename = os.path.dirname(file_other) + "/mismatches_vs_GR_v3.dat"
    header = "          Mass         Mismatch          <h1,h2>          <h1,h1>          <h2,h2>  Best Time Shift   Peak Time Diff"
    np.savetxt(filename, mismatches, header=header, fmt="%.10e")


if __name__ == "__main__":

    filename = "strain_22_cutoff_low_0.015_at_infinity.dat"
    # folder_source = "g3_high/run_g3_negative_v2_h112_higher_value/"
    folder_source = "GR/run_GR_v2_h112/"

    folders_compare = [
        "g2_0.02/run_g2_v2_h112/",
        "g3_high/run_g3_negative_v2_h112_higher_value/",
    ]

    for folder in folders_compare:
        compute_all_mismatch_data(folder_source + filename, folder + filename)


# nuances that may change the mismatch:
# - taper window percentage
# - junk radiation cutoff
# - final end time of the wave
# - solar masses
# - cut off high frequency
# - cut off low frequency
# - discretization of the time series (interpolation with more points)
