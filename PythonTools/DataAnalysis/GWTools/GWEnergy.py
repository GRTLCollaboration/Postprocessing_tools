
"""
Compute GWEnergy from Weyl4ExtractionOut files
"""

import sys
sys.path.append("../Base/")

import SmallDataIOReader, IntegrationMethod, DataIntegration
import math

file_prefix = "run_g2_v3_h128/weyl4/Weyl4ExtractionOut_*.dat"
reader = SmallDataIOReader.FileSet(file_prefix, read = False)

# # General info about files
# file0 = reader.getFile(0)
# file0.read()
# print(file0.getName())
# print(file0.wasRead())
# print(file0.getHeaders())
# print(file0.numBlocks())
# print(file0.numLabels())
# print(file0.numValues())
# print(file0.numRows())
# print(file0.numColumns())
# print(file0.getData())

# Integrate all files over time (for each point in space)
method = IntegrationMethod.midpoint
weyl4_time_integrated = DataIntegration.integrate_in_time(reader, method, verbose = True, max_steps = -1)

# add column with abs(Weyl4)
weyl4_time_integrated.addColumn(lambda row, map: row[map['Weyl4_Re'][1]] ** 2 + row[map['Weyl4_Im'][1]] ** 2)
# ignore other columns
weyl4_time_integrated.removeColumn(2, 3)

# Now integrate the resulting time integration in space
weyl4_time_and_space_integrated = DataIntegration.integrate_in_space_2d(weyl4_time_integrated,
                                        [IntegrationMethod.simpson, IntegrationMethod.trapezium],
                                        [False, True],
                                        'r',
                                        DataIntegration.spherical_area_element)

# don't forget those pi's
for result in weyl4_time_and_space_integrated:
    result[1][0] /= (16. * math.pi)

print(weyl4_time_and_space_integrated)

# Quick attempt just with 22 mode:

# filename = "data/Weyl_integral_22.dat"
# mode22 = SmallDataIOReader.File(filename)

# accumulated = integrate_in_time(mode22[0], method, accumulate = True)
# print(accumulated)
# print(integrate_in_time(mode22[0], method))
# print([0][:2])
