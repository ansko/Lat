#!/usr/bin/env python3

import copy
import math
import Lat.Utilities as uti


#######################################################
#
# reads dump file, return list stress_dump[stress_x][z]
#

def read_dump(datafile, mode, stress_dump):
    f=open(datafile, 'r')
    i=0

    if mode not in [2, 7]:
        print("Incorrect mode")
        return None

    for string in f:
        i += 1
        if i < 10:
            continue
        string_splitted=string.split()
        (stress_dump[int(string_splitted[0])])[0] = float(string_splitted[1])
        (stress_dump[int(string_splitted[0])])[1] = float(string_splitted[2])

    return stress_dump


#####################################
#
# converts atomic stresses to pressure
# for comp + 2 phases
#

def pressure_two_phases(stress_dump):
    cellvolume = 310321
    softvolume = 243321

    pres_comp = pres_mmt = pres_soft = 0

    for i in range(1, 31321):
        pres_comp += stress_dump[i][1]
        if (i < 721 or
           (i > 3480 and i < 4201) or
           (i > 6960 and i < 7681) or 
           (i > 10440 and i < 11160) or
           (i > 13920 and i < 14641) or
           (i > 17420 and i < 18121) or
           (i > 20880 and i < 21601) or
           (i > 24360 and i < 25081) or (i > 27840 and i < 28561)):
            pres_mmt += stress_dump[i][1]
        else:
            pres_soft += stress_dump[i][1]

    return [pres_comp / cellvolume, 
            pres_mmt / (cellvolume - softvolume),
            pres_soft / softvolume]


#######################################
#
# converts atomic stresses to pressure
# for for 7 layers
#

def pressure_layers(stress_dump):
    cellvolume = 310321
    softvolume = 243321
#    layers_bounds=[60, 64.5, 69.5, 73.5, 78, 82.5, 87, 92]
    layers_bounds=[13, 18, 22.5, 27, 32, 36, 40.5, 45]

    pres_layers = [0, 0, 0, 0, 0, 0, 0, 0]

    for i in range(1, 31321):
        if layers_bounds[0] < stress_dump[i][0] < layers_bounds[1]:
            pres_layers[0] += stress_dump[i][1]
        elif layers_bounds[1] < stress_dump[i][0] < layers_bounds[2]:
            pres_layers[1] += stress_dump[i][1]
        elif layers_bounds[2] < stress_dump[i][0] < layers_bounds[3]:
            pres_layers[2] += stress_dump[i][1]
        elif layers_bounds[3] < stress_dump[i][0] < layers_bounds[4]:
            pres_layers[3] += stress_dump[i][1]
        elif layers_bounds[4] < stress_dump[i][0] < layers_bounds[5]:
            pres_layers[4] += stress_dump[i][1]
        elif layers_bounds[5] < stress_dump[i][0] < layers_bounds[6]:
            pres_layers[5] += stress_dump[i][1]
        elif layers_bounds[6] < stress_dump[i][0] < layers_bounds[7]:
            pres_layers[6] += stress_dump[i][1]
        else:
            pres_layers[7] += stress_dump[i][1]

    return [pres_layers[0] / cellvolume *
            (layers_bounds[1] - layers_bounds[0]) / 41, 

            pres_layers[1] / cellvolume *
            (layers_bounds[2] - layers_bounds[1]) / 41,

            pres_layers[2] / cellvolume *
            (layers_bounds[3] - layers_bounds[2]) / 41,

            pres_layers[3] / cellvolume *
            (layers_bounds[4] - layers_bounds[3]) / 41,

            pres_layers[4] / cellvolume *
            (layers_bounds[5] - layers_bounds[4]) / 41,

            pres_layers[5] / cellvolume *
            (layers_bounds[6] - layers_bounds[5]) / 41,

            pres_layers[6] / cellvolume *
            (layers_bounds[7] - layers_bounds[6]) / 41,

            pres_layers[7]]


#############################################
#
# approximates by first term of Fourier seria
#

def appro_fourier(pressure, period):
    meanings = [[]]
    i = a0 = a1 = a2 = b1 = b2 = 0
    suma0 = suma1 = suma2 = sumb1 = sumb2 = magnitude = 0

    for i in range(1, period):
        meanings.append(float(pressure[i][0])) #0 - comp, 1 - mmt, 2 - soft
                                               #or for layers
        suma0 += float(meanings[i])
        suma1 += float(meanings[i]) * math.cos(2 * math.pi / period * (i + 1))
        sumb1 += float(meanings[i]) * math.sin(2 * math.pi / period * (i + 1))
        suma2 += float(meanings[i]) * math.cos(4 * math.pi / period * (i + 1))
        sumb2 += float(meanings[i]) * math.sin(4 * math.pi / period * (i + 1))

    a0 = suma0 / period
    a1 = 2 * suma1 / period
    a2 = 2 * suma2 / period
    b1 = 2 * sumb1 / period
    b2 = 2 * sumb2 / period
    magnitude = 5 * math.sqrt(a1**2 + b1**2)
   # print(magnitude)
    return magnitude
   # return None
