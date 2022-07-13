
"""
Convert Weyl4 to Strain (complex values) using FFT

NB: Time domain direct integral has a lot of error noise when I tried
(check Section F of https://arxiv.org/abs/2112.15529)
"""

import sys, os, glob
sys.path.append("../Base/")
sys.path.append("../")
from utils import *
import SmallDataIOReader

import numpy as np
import matplotlib.pyplot as plt

# cut offs are normalized by 2pi/(N*dt), so they vary between 1/N and 1/2
# where N is the length of the zero-padded data
# useful: our Comparison paper
# useful: arxiv 1006.1632
# useful: my paper (Cubic Horndeski Binaries)
def weyl_to_strain(filepath, folder_out = None, name_prefix = "strain", output = None, cut_off_low = -1, cut_off_high = 1e10, junk_time=0):
    assert "Weyl4_mode" in os.path.basename(filepath)
    assert cut_off_low > 0

    mode = SmallDataIOReader.File(filepath).getBlock(0)
    mode_data = np.transpose(mode.getData())

    # chop off junk radiation
    mode_data = mode_data[:, mode_data[0]>=junk_time ]

    times = mode_data[0]
    dt = times[1] - times[0]

    freq_w_low = 2. * math.pi * cut_off_low
    freq_w_high = 2. * math.pi * cut_off_high

    strain_data = [times]

    file_headers = mode.getHeaders()
    # original files vs extrapolated to infinity files
    final = len(mode_data) if len(file_headers)>1 else (len(mode_data)+1)//2

    for i in range(1, final, 2):
        # time-domain to FFT
        freqs, mode_FFT = FFT(mode_data[i] + 1j * mode_data[i+1])

        # the time-domain integral is just a division in Fourier space
        freqs_w = 2. * math.pi * freqs / dt
        mode_FFT_strain = - mode_FFT / freqs_w**2

        if i==1:
            print("Lowest frequency ratio:", freqs[0], f"(angular frequency {freqs_w[0]})")
            print("Using low cutoff:", cut_off_low, f"(angular frequency {freq_w_low})")

        mode_FFT_strain[np.abs(freqs_w) < freq_w_low] = - mode_FFT[np.abs(freqs_w) < freq_w_low] / freq_w_low**2
        mode_FFT_strain[np.abs(freqs_w) > freq_w_high] = - mode_FFT[np.abs(freqs_w) > freq_w_high] / freq_w_high**2

        # FFT back to time-domain
        mode_strain = iFFT(mode_FFT_strain)

        # plus strain is minus the real part
        mode_strain = -np.conjugate(mode_strain)
        strain_data.append(np.real(mode_strain[:len(times)]))
        strain_data.append(np.imag(mode_strain[:len(times)]))

    # original files vs extrapolated to infinity files:
    if len(file_headers) > 1:
        headers_1 = '\t\t\t\t'.join(mode.getHeaders()[0]).replace("integral Re", "h_+").replace("integral Im", "h_x")
        headers_2 = '\t\t  '.join(mode.getHeaders()[1])
        header = headers_1 + "\n" + headers_2
    else:
        header = '\t\t\t\t'.join(mode.getHeaders()[0]).replace("Re(Weyl)", "h_+").replace("Im(Weyl)", "h_x")

    folder_out = folder_out if folder_out is not None else (os.path.dirname(filepath) + "/")
    output = folder_out + os.path.splitext(os.path.basename(filepath))[0].replace("Weyl4_mode", "strain") + f"_cutoff_low_{cut_off_low}.dat"
    np.savetxt(output, np.transpose(strain_data), header = header, fmt="%.10e")


if __name__ == "__main__":
    # important - low frequency cut off (play with it until you see you minimize noise without losing data)
    cut_off_low = 0.015
    # start time, after junk radiation
    junk_time = 0
    
    # weyl_to_strain all modes on these folders:
    folders = ["./g1/"]

    for folder in folders:
        print("ON FOLDER:", folder)

        file_list = glob.glob(f"{folder}Weyl4_mode_??.dat")

        for file_path in file_list:
            print("On file:", file_path)
            weyl_to_strain(file_path, folder, cut_off_low=cut_off_low, junk_time=junk_time)

