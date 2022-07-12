
import sys
sys.path.append("../Base/")
sys.path.append("../")
from utils import *
import SmallDataIOReader

from scipy.interpolate import make_interp_spline
import numpy as np

# get the time of merger for the following files:
# 2nd parameter is the spacetime mass (it helps aligning the waves)
files = [
    ("./run_g3_positive_v2_h112/Weyl_integral_22.dat", initial_mass_g3_0_01_v2),
    ("./run_g3_negative_v2_h112/Weyl_integral_22.dat", initial_mass_g3_m0_01_v2),
]

##########

def get_t_merger(file_path, mass):

    file = SmallDataIOReader.File(file_path)

    file.complexifyColumns()
    data = np.array(file.getData())
    time = np.real(data[:,0])
    amps = np.transpose(np.abs(data[:,1:]))

    radii = get_radii(file_path, verbose=False)
    tortoise_rs = [tortoise_radius(radius, mass) for radius in radii]

    weyls_interp = [make_interp_spline(time, amp, k=3) for amp in amps]
    times_new = np.linspace(time[0], time[-1], 1000000)
    weyls_full = [interp(times_new) for interp in weyls_interp]

    max_time = [times_new[np.argmax(weyl)] for weyl in weyls_full]

    time_retarded = [time - r for time, r in zip(max_time, tortoise_rs)]

    mean = np.mean(time_retarded)
    error = np.max(np.abs(time_retarded-mean))
    
    print(file_path)
    # print(mean, "+-", error)
    print(f"{mean:.2f} +- {error + 0.005:.2f}") # round error up to 2nd decimal place

    return mean, error

for file in files:
    get_t_merger(file[0], file[1])
