
"""
Use this file to integrate data in time or space
Useful functions are:
- 'integrate_in_time'
- 'integrate_in_space_2d'
They can integrate a list/array, a SmallDataIO Block, File or set of Files

E.g. To integrate over the AH:
DataIntegration.integrate_in_space_2d(coordsAHFile,
                                    [IntegrationMethod.trapezium, IntegrationMethod.trapezium],
                                    [False, True],
                                    'r',
                                    DataIntegration.spherical_area_element)

E.g. To integrate Weyl Extraction files over time (integrate a set of files, for all spatial points)
weyl4_time_integrated = DataIntegration.integrate_in_time(fileSet, someMethod, verbose = True, max_steps = -1)

"""

import SmallDataIOReader, IntegrationMethod
import copy # use deepcopy
import numpy
import math

def computeDt(data, x):
    """ x only needed if data is a ndarray or list"""
    if isinstance(data, SmallDataIOReader.Block):
        return data[1][0] - data[0][0]
    elif isinstance(data, SmallDataIOReader.FileSet):
        if not data[0].wasRead(): data[0].read()
        if not data[1].wasRead(): data[1].read()
        t1 = data[1][0].getHeaderValue('time')
        assert(t1 != None)
        return t1 - data[0][0].getHeaderValue('time')
    elif isinstance(data, SmallDataIOReader.File):
        assert("Are you sure you want to integrate a File and not a Block?" == "")
    elif isinstance(data, numpy.ndarray) or isinstance(data, list):
        return x[1] - x[0]
    else:
        assert("Not implemented yet" == "")

def initialState(data):
    if isinstance(data, SmallDataIOReader.Block):
        if isinstance(data[0], numpy.ndarray):
            return data[0][data.numLabels():]*0 # return only the labels
        elif isinstance(data[0], list):
            return numpy.array(data[0][data.numLabels():])*0 # return only the labels
    elif isinstance(data, SmallDataIOReader.FileSet):
        return data[0]*0 # return a zero-ed file with labels preserved
    elif isinstance(data, SmallDataIOReader.File):
        assert("Are you sure you want to integrate a File and not a Block?" == "")
    elif isinstance(data, numpy.ndarray):
        return data[0]*0 # assume no labels
    elif isinstance(data, list):
        return numpy.array(data[0])*0 # assume no labels
    else:
        assert("Not implemented yet" == "")

def prepareForIntegral(data_part, data):
    if isinstance(data_part, numpy.ndarray): # part of a Block
        if isinstance(data, numpy.ndarray):
            return data_part # all good, assume no labels
        else:
            return data_part[data.numLabels():] # line without labels
    elif isinstance(data_part, list): # part of a Block that was not numpy'ed
        if isinstance(data, list):
            return numpy.array(data_part); # all good, assume no labels
        else:
            return numpy.array(data_part[data.numLabels():]) # line without labels
    elif isinstance(data_part, SmallDataIOReader.File): # part of a FileSet
        if not data_part.wasRead(): data_part.read()
        return data_part
    elif isinstance(data_part, SmallDataIOReader.Block): # part of a File
        assert("Are you sure you want to integrate a File and not a Block?" == "")
    else:
        assert("Not implemented yet" == "")

# 'accumulate' flags only work with Trapezium or Midpoint method, as simpson can't accumulate in the middle of its 2 weight step
def integrate_in_time(data, method, verbose = False, max_steps = None, accumulate = False, accumulate_twice = False, x = None):
    # checks
    assert accumulate or not accumulate_twice # can't have accumulate = False and accumulate_twice = True
    assert isinstance(method, IntegrationMethod.IntegrationMethod)
    assert (isinstance(data, SmallDataIOReader.Block) or
           isinstance(data, SmallDataIOReader.File) or
           isinstance(data, SmallDataIOReader.FileSet) or
           isinstance(data, numpy.ndarray) or isinstance(data, list))

    is_periodic = False
    dt = computeDt(data, x)
    integral = initialState(data)
    times = [0]
    accumulated_integral = [integral]

    # for accumulate_twice:
    double_integral = initialState(data)
    accumulate_twice_integral = [integral]

    if max_steps is not None:
        data_cut = data[:max_steps]
    else:
        data_cut = data
    total_steps = len(data_cut)

    assert(method.isValid(total_steps, is_periodic))

    for step, part in enumerate(data_cut):
        if verbose:
            print("Step = %d / %d" % (step+1, total_steps))
        weight = method.weight(step, total_steps, is_periodic)
        to_add = prepareForIntegral(part, data)
        integral += dt * weight * to_add
        if accumulate:
            accumulated_integral.append(copy.deepcopy(integral))
            times.append(times[-1]+dt)
            if accumulate_twice:
                double_integral += dt * weight * integral
                accumulate_twice_integral.append(copy.deepcopy(double_integral))

    if accumulate_twice:
        return times, accumulated_integral, accumulate_twice_integral
    if accumulate:
        return times, accumulated_integral
    return integral

def get_2d_surface_dimensions(block):
    # checks
    assert(isinstance(block, SmallDataIOReader.Block))
    assert(block.numLabels() in (2, 3)) # we are in 2d files (coords files have radius, so 3 labels)

    data = block.getData()
    n = len(data)
    assert(n >= 2)

    index_u = 1
    while index_u < len(data) and data[index_u][0] == data[0][0]:
        index_u = index_u + 1
    assert(index_u < len(data))
    assert(n // index_u == n / index_u)
    du = data[index_u][0] - data[0][0]

    index_v = 1
    while index_v < len(data) and data[index_v][1] == data[0][1]:
        index_v = index_v + 1
    assert(index_v < len(data))
    assert(n // index_v == n / index_v)
    dv = data[index_v][1] - data[0][1]

    if index_u > index_v:
        return (n // index_u, index_u, du, dv)
    else:
        return (index_v, n // index_v, du, dv)

def unit_area_element(surface_param_value, u, v):
    return 1.
def spherical_area_element(surface_param_value, u, v):
    return (surface_param_value ** 2) * math.sin(u)

def integrate_in_space_2d(data, methods, is_periodic, header_label, area_element):
    # checks
    assert(len(methods) == 2 and len(is_periodic) == 2)
    assert(isinstance(methods[0], IntegrationMethod.IntegrationMethod))
    assert(isinstance(methods[1], IntegrationMethod.IntegrationMethod))

    if isinstance(data, SmallDataIOReader.Block):
        surface_param_value = data.getHeaderValue(header_label)
        surface_param_value_found = surface_param_value

        (num_u, num_v, du, dv) = get_2d_surface_dimensions(data)
        integral = initialState(data)

        assert(methods[0].isValid(num_u, is_periodic[0]))
        assert(methods[1].isValid(num_v, is_periodic[1]))

        for iu in range(num_u):
            u = data[iu * num_v][0]
            inner_integral = initialState(data)
            for iv in range(num_v):
                v = data[iv][1]

                if surface_param_value == None:
                    index = data.findHeader(header_label)
                    assert index is not None
                    surface_param_value_found = data[iv][index[1]]

                part = data[iu * num_v + iv]
                weight_v = methods[1].weight(iv, num_v, is_periodic[1])
                to_add = prepareForIntegral(part, data)
                inner_integral += dv * weight_v * to_add * area_element(surface_param_value_found, u, v)
                if numpy.isnan(inner_integral).any():
                    print(f"Nans found for u={u}, v={v}. Setting them to 0")
                    inner_integral[numpy.isnan(inner_integral)] = 0
                if (numpy.abs(inner_integral) > 10*10).any():
                    print(f"Big values found for u={u}, v={v}. Setting them to 0")
                    inner_integral[numpy.abs(inner_integral) > 10*10] = 0

            weight_u = methods[0].weight(iu, num_u, is_periodic[0])
            integral += du * weight_u * inner_integral

        return (surface_param_value, integral) if surface_param_value is not None else integral
    elif isinstance(data, SmallDataIOReader.File):
        return numpy.array([integrate_in_space_2d(block, methods, is_periodic, header_label, area_element) for block in data.blocks])
    elif isinstance(data, SmallDataIOReader.FileSet):
        return numpy.array([integrate_in_space_2d(file, methods, is_periodic, header_label, area_element) for file in data.files])
    else:
        raise NotImplementedError("Type " + str(type(data)) + " not implemented yet.")
