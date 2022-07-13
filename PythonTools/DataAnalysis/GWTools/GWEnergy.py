
"""
Compute GWEnergy from Weyl4ExtractionOut files

NB: following 8.9.32 of Alcubierre

"""

import sys, glob
sys.path.append("../Base/")
sys.path.append("../")
from utils import *
import SmallDataIOReader, IntegrationMethod, DataIntegration

import math
import numpy as np
import multiprocessing

method = IntegrationMethod.midpoint

file_prefix = base_folder + "Weyl4ExtractionOut_*.dat"
skip = 1

reader = SmallDataIOReader.FileSet(file_prefix, read = False, skip=skip)

# read files in parallel
numFiles = reader.numFiles()
def readFileAndCutBlock(obj):
    # unpack
    i, file = obj
    print(f"Reading {i+1}/{numFiles}")
    file.read()
    # only keep r = 50, r = 60
    # file.removeBlock(2, file.numBlocks()-1)
    return file

pool_obj = multiprocessing.Pool()
answer = pool_obj.map(readFileAndCutBlock, enumerate(reader))
pool_obj.close()
pool_obj.join()

reader = SmallDataIOReader.FileSet(files = answer)

times, weyl4_time_integrated_all_time = DataIntegration.integrate_in_time(reader, method, verbose = True, accumulate = True)

def spatial_integration(obj):
    # unpack
    i, weyl4_time_integrated = obj

    print(f"Spacial integration {i+1}/{numFiles}")
    weyl4_time_integrated.addColumn(lambda row, map: row[map['Weyl4_Re'][1]] ** 2 + row[map['Weyl4_Im'][1]] ** 2)
    weyl4_time_integrated.removeColumn(2, 3)

    weyl4_time_and_space_integrated = DataIntegration.integrate_in_space_2d(weyl4_time_integrated,
                                            [IntegrationMethod.simpson, IntegrationMethod.trapezium],
                                            [False, True],
                                            'r',
                                            DataIntegration.spherical_area_element)

    power = []
    for result in weyl4_time_and_space_integrated:
        power.append( result[1][0] / (16. * math.pi) )

    return power

pool_obj = multiprocessing.Pool()
power_all_time = pool_obj.map(spatial_integration, enumerate(weyl4_time_integrated_all_time))
pool_obj.close()
pool_obj.join()

print("DONE COMPUTING POWER")

# print(power_all_time)

times, energies = DataIntegration.integrate_in_time(power_all_time, method, x=times, verbose=True, accumulate=True)

radii = [str(block.getHeaderValue('r')) for block in reader[0]]
header = 'time\t\tr = ' + '\t\tr = '.join(radii)

output = np.transpose(np.concatenate([[times], np.transpose(energies)], axis=0))

np.savetxt(f"GW_energies.dat", output, header=header, fmt="%.10f")

# mine for r=50 and r=60:
# 0.11218792
# 0.06091224
# Miren C++:
# 0.11209019876
# 0.060865933070
# error: <0.1%
