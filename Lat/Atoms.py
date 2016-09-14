#!/usr/bin/env python3


import math
import copy
import Lat.Utilities as uti


#
# reads atoms from .data, return list with atoms:
# [atomtype, atom_x, atom_y, atom_z]
#

def read_atoms(datafile):
    f = open(datafile, 'r')
    atoms=[[]]
    start=0

    for i in range(1, 31321):
        atoms.append([0, 0, 0, 0])

    i = 0

    for line in f:
        i += 1
        if line.startswith("Atoms"):	#start reading here
            start = i + 1
        if line.startswith("Velocities"): # finsih reading here
            break
        if i > start > 0:
            line_splitted = line.split()
            if len(line_splitted) == 0:
                continue
            atoms[int(line_splitted[0])] = [float(line_splitted[2]),
                                            float(line_splitted[4]),
                                            float(line_splitted[5]),
                                            float(line_splitted[6])]

    return atoms


################################################
#
# calculates atomic displacement in different phases
#

def phase_mobility(atoms_old, atoms_new, disp_vector):
    displacement=[0, 0, 0]
    bounds=[3.6574177352488846e-01, 9.3767139801735425e+01,
            -1.7173796716130170e+00, 7.8533787176849529e+01,
            4.2580775568494040e+00, 4.6615800720898385e+01]

    for i in range(1, len(atoms_old)):
        delta=atoms_new[i][1] - atoms_old[i][1]
        delta2=bounds[1] - bounds[0] - abs(delta)
        if abs(delta) <= abs(delta2):
            disp_vector[i-1][0] += delta
        else:
            disp_vector[i-1][0] += math.copysign(delta2, delta)

        delta=atoms_new[i][2] - atoms_old[i][2]
        delta2=bounds[3] - bounds[2] - abs(delta)
        if abs(delta) <= abs(delta2):
            disp_vector[i-1][1] += delta
        else:
            disp_vector[i-1][1] += math.copysign(delta2, delta)

        delta=atoms_new[i][3] - atoms_old[i][3]
        delta2=bounds[5] - bounds[4] - abs(delta)
        if abs(delta) <= abs(delta2):
            disp_vector[i-1][2] += delta
        else:
            disp_vector[i-1][2] += math.copysign(delta2, delta)

    return disp_vector


################################################
#
# calculates and prints atomic mobility in different phases
#

def print_phase_mobility(disp_vector, mode):

    atoms_disp=atoms_disp_x=atoms_disp_y=atoms_disp_z=[0, 0, 0]

    for i in range(1, 31320):
        atoms_disp_x[uti.define_phase(i) - 1] += disp_vector[i-1][0]**2
        atoms_disp_y[uti.define_phase(i) - 1] += disp_vector[i-1][1]**2
        atoms_disp_z[uti.define_phase(i) - 1] += disp_vector[i-1][2]**2

    if mode==0:
        for coordinate in (atoms_disp_x, atoms_disp_y, atoms_disp_z): #(<x^2> (t))
            coordinate=[coordinate[0] / 720,
                        coordinate[1] / 840,
                        coordinate[2] / 1920]
    elif mode==1:
        for coordinate in (atoms_disp_x, atoms_disp_y, atoms_disp_z): #(<x>^2 (t))
            coordinate=[coordinate[0] / 720**2,
                        coordinate[1] / 840**2,
                        coordinate[2] / 1920**2]

    atoms_disp[0]=(atoms_disp_x[0] + atoms_disp_y[0] + atoms_disp_z[0])
    atoms_disp[1]=(atoms_disp_x[1] + atoms_disp_y[1] + atoms_disp_z[1])
    atoms_disp[2]=(atoms_disp_x[2] + atoms_disp_y[2] + atoms_disp_z[2])

    print(atoms_disp[0], atoms_disp[1], atoms_disp[2])


####################
#
# prints atoms array
#

def print_atoms(atoms):
    for i in range(1, 31321):
        print(atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3])

    return None

############################
#
# builds 4 density profiles:
# mmt
# modifier
# polymer
# whole composite
#

def components_density_profile(atoms, distribution, scale):
   fluctuations = 100 # to handle with drift to negative idicies

   for i in range(1, len(atoms)):
      distribution[fluctuations + int(atoms[i][3] * scale)][0] += 1
      distribution[fluctuations + 
                   int(atoms[i][3] * scale)][uti.define_phase(i)] += 1

   return distribution


###############################################
#
# quadric estimation of surfactant distribution
#

def appro_modifier(distribution):
    estimation = 0
    surfactant=[]

    for i in range(len(distribution)):
        if distribution[i][2] != 0:
            surfactant.append(distribution[i][2])

    surfactant_rev=copy.deepcopy(surfactant)
    surfactant_rev.reverse()

    for i in range(int(len(surfactant)/2) + 1):
        surfactant[i] += surfactant_rev[i]

    for i in range(int(len(surfactant)/2) + 1):
        estimation += surfactant[i]**2

    return estimation
