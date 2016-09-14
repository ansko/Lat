#!/usr/bin/env python3

#################################
#
# main functions are in this file
#

import copy
import math
import threading
import Lat.Utilities as uti
import Lat.Bonds as Bond
import Lat.Atoms as Atom
import Lat.Stress as Stress

############################################
#
# calculates curve (<x^2> (t)) or (<x>^2(t))
#

def atoms_displacement():
    datafile=""
    bondsfile="bonds_seg"
    j=1

    foldernames=[]
    foldernames2=[]
    filenames=[]

    foldernames2.append("0 - relaxation (1410408)")

    foldernames.append("0 - relaxation (1410408)")
    foldernames.append("1st wiggle cycle (1414047)")
    foldernames.append("2nd wiggle cycle (1426324)")
    foldernames.append("3rd wiggle cycle corrupted (1427112)")
    foldernames.append("4th wiggle cycle after corrupted (1427697)")
    foldername="3rd wiggle cycle corrupted (1427112)"

    disp_vector=[[0, 0, 0] for i in range(31320)] 
                                        # add displacements every step

    k=0

    for foldername in foldernames:
        for j in range(1, 51):
            filenames.append("".join(("seg/", foldername, "/co.",
                                       str((j)*5), "0000.data")))

    for j in range(0, 249):
        datafile_old=filenames[j]
        datafile_new=filenames[j+1]

        atoms_old=atoms.read_Atom(datafile_old)
        atoms_new=atoms.read_Atom(datafile_new)

        disp_vector=Atom.phase_mobility(atoms_old, atoms_new,
                                        disp_vector)
        print(k, end = " ")
        Atom.print_phase_mobility(disp_vector, 0)
                                            #mode=0 for <x^2>, 1 for <x>^2

        datafile_old=""
        datafile_new=""
        k+=1

    return None

################################
#
# calculates bond or bond 1-3 orientation parameter
#

def bond_orientation_parameter():
#    bondsfile="bonds_seg" #bonds
    bondsfile="angles_seg" #bonds 1-3
    atomsfile="seg/0 - relaxation (1410408)/co.2500000.data"

    bonds=Bond.read_bonds(bondsfile)
    orientation=None

    foldernames=[]
    foldernames.append("0 - relaxation (1410408)")
    foldernames.append("1st wiggle cycle (1414047)")
    foldernames.append("2nd wiggle cycle (1426324)")
    foldernames.append("3rd wiggle cycle corrupted (1427112)")
    foldernames.append("4th wiggle cycle after corrupted (1427697)")

    k = 0

    for foldername in foldernames:
        for j in range(1, 51):
            k += 1
            atomsfile="".join(("seg/", foldername, "/co.",
                                str(j*5), "0000.data"))
            atoms=Atom.read_atoms(atomsfile)
            orientation=Bond.bond_orientation(atoms, bonds,
                                              k, orientation, 1, 5)

    for i in range(1, len(orientation)):
        if orientation[i][1] != 0:
            print(i, 200 * orientation[i][0] / orientation[i][1], 
                  orientation[i][1] / (50 * len(foldernames)))

    return None

############################
#
# To get phase distribution
#

def phase_distribution():
    stepsnumber=50
    scale=3


    foldernames=[]
    foldernames.append("0 - relaxation (1410408)")
    foldernames.append("1st wiggle cycle (1414047)")
    foldernames.append("2nd wiggle cycle (1426324)")
    foldernames.append("3rd wiggle cycle corrupted (1427112)")
    foldernames.append("4th wiggle cycle after corrupted (1427697)")

    distribution=[[0, 0, 0, 0] for i in range(1, 31320*scale)]
                                                 #composite, mmt, mod, poly
    for foldername in foldernames:
        for i in range(1, stepsnumber):
            atoms=Atom.read_atoms("".join(("seg/",
                                           foldername,
                                           "/co.",
                                           str(i*5),
                                           "0000.data")))
            distribution=Atom.components_density_profile(atoms,
                                                         distribution,
                                                         scale)

    for i in range(1, len(distribution)):
        if(distribution[i][0] != 0):
            print((i - 100) / 3, distribution[i][0] / stepsnumber,
                     distribution[i][1] / stepsnumber,
                     distribution[i][2] / stepsnumber,
                     distribution[i][3] / stepsnumber)

    return None


#####################################
#
# To estimate surfactant distribution
#

def surfactant_distribution():
    stepsnumber = 50
    scale = 3


    foldernames_seg = []  
#    foldernames_seg.append("-6 - 1332494 - 2.5ns")
    foldernames_seg.append("-5 - 1338014 - 2.5ns more")
    foldernames_seg.append("-4 - 1343289 - 2.5ns more")
    foldernames_seg.append("-3 - 1352547 - 2.5 ns more")
    foldernames_seg.append("-2 - 1388142 - 2.5 ns more")
    foldernames_seg.append("-1 - 1404644")
    foldernames_seg.append("0 - relaxation (1410408)")
    foldernames_seg.append("1st wiggle cycle (1414047)")
    foldernames_seg.append("2nd wiggle cycle (1426324)")
    foldernames_seg.append("3rd wiggle cycle corrupted (1427112)")
    foldernames_seg.append("4th wiggle cycle after corrupted (1427697)")
    foldernames = copy.deepcopy(foldernames_seg)

    distr = [[0, 0, 0, 0] for i in range(1, 31320 * scale)]
                                            #composite, mmt, mod, poly

    f = open('dist_seg', 'w')

    for foldername in foldernames:
        for i in range(1, 51, 5):
            distribution=copy.deepcopy(distr)
            atoms=Atom.read_atoms("".join(("seg/", foldername,
                                   "/co.", str(i * 5), "0000.data")))
            distribution = Atom.components_density_profile(atoms,
                                                           distribution,
                                                           scale)

            estimation = Atom.appro_modifier(distribution)

            f.write(str(estimation / 1000000) + '\n')

    foldernames_mix = []
    foldernames_mix.append("-5 - 1388141")
    foldernames_mix.append("-4 - 1393565")
    foldernames_mix.append("-3 - 1397381")
    foldernames_mix.append("-2 - 1400005")
    foldernames_mix.append("-1 - 1404643")
    foldernames_mix.append("0 - relaxation (1410378)")
    foldernames_mix.append("1st wiggle cycle (1414048 )")
    foldernames_mix.append("2nd wiggle cycle (1426323)")
    foldernames_mix.append("3rd wiggle cycle(1427113)")
    foldernames_mix.append("4th wiggle cycle (1427696)")
    foldernames = copy.deepcopy(foldernames_mix)

    f.close()

    f = open('dist_mix', 'w')

    for foldername in foldernames:
        for i in range(1, 51, 5):
            distribution = copy.deepcopy(distr)
            atoms = Atom.read_atoms("".join(("mix/", foldername,
                                             "/co.", str(i * 5),
                                             "0000.data")))
            distribution = Atom.components_density_profile(atoms,
                                                           distribution,
                                                           scale)

            estimation = Atom.appro_modifier(distribution)

            f.write(str(estimation / 1000000) + '\n')

    f.close()

    return None

###########################################
#
# reads dumps and makes stress-strain curve
#

def stress_strain_two_phases():
    foldernames_mix = []
    foldernames_mix.append("1st wiggle cycle (1414048 )")
    step = 100

    pressure = []

    stress_dump = [[0, 0] for i in range(31321)]

    for foldername in foldernames_mix:
        for i in range(0, 25001, step):
            datafile = "".join(("mix/", foldername, "/Dumps/ALLstress.",
                               str(i * 100)))
            stress_dump = Stress.read_dump(datafile, 2, stress_dump)
#            pressure.append(Stress.pressure_two_phases(stress_dump))
            pressure.append(Stress.pressure_layers(stress_dump))

    f = open('pressure', 'w')
    for i in range(len(pressure)):
        f.write(str(i) + " " +
                str(pressure[i][0]) + " " +#phases
                str(pressure[i][1]) + " " +#
                str(pressure[i][2]) + " " +##layers
                str(pressure[i][3]) + " " +##
                str(pressure[i][4]) + " " +##
                str(pressure[i][5]) + " " +##
                str(pressure[i][6]) + " " +##
                str(pressure[i][7]) + "\n")##
#                str(pressure[i][2]) + "\n")#
    f.close()

    magnitude = Stress.appro_fourier(pressure, int(25001/step)-1)

    print(magnitude)

    return None


###########################################
#
# reads dumps and makes stress-strain curve
# (with threads)
#

def stress_strain_two_phases_p(n, datafile, step):
    foldernames_mix = []
    foldernames_mix.append("1st wiggle cycle (1414048 )")

    pressure0 = []
    pressure1 = []
    pressure2 = []
    pressure3 = []

    stress_dump = [[0, 0] for i in range(31321)]

    for foldername in foldernames_mix:
        if n==0:
            for i in range(0, 6251, step):
                pressure0 = for_stress_strain_two_phases_p(i, foldername,
                                                          stress_dump,
                                                          pressure0)
            for_stress_strain_two_phases_p2(datafile, 0, 6251, 
                                            pressure0, step)
        if n==1:
            for i in range(6250, 12501, step):
                pressure1 = for_stress_strain_two_phases_p(i, foldername,
                                                          stress_dump,
                                                          pressure1)
            for_stress_strain_two_phases_p2(datafile, 6250, 12501, 
                                            pressure1, step)
        if n==2:
            for i in range(12500, 18751, step):
                pressure2 = for_stress_strain_two_phases_p(i, foldername,
                                                          stress_dump,
                                                          pressure2)
            for_stress_strain_two_phases_p2(datafile, 12500, 18751,
                                            pressure2, step)
        if n==3:
            for i in range(18750, 25001, step):
                pressure3 = for_stress_strain_two_phases_p(i, foldername,
                                                          stress_dump,
                                                          pressure3)
            for_stress_strain_two_phases_p2(datafile, 18750, 25001,
                                            pressure3, step)


####################################################
#
# accessory functions for stress_strain_two_phases_p
#

def for_stress_strain_two_phases_p(i, foldername, stress_dump, pressure):
    datafile = "".join(("mix/", foldername, "/Dumps/ALLstress.",
                         str(i * 100)))
    stress_dump = Stress.read_dump(datafile, 2, stress_dump)
    pressure.append(Stress.pressure_two_phases(stress_dump))

    return pressure


def for_stress_strain_two_phases_p2(datafile, start, stop, pressure, step):
    f = open(datafile, 'w')
    for i in range(len(pressure)):
        f.write(str(start//step + i) + " " +
                str(pressure[i][0]) + " " +#phases
                str(pressure[i][1]) + " " +#
                str(pressure[i][2]) + "\n")#
    f.close()

    return None


####################################
#
# to call stress_strain_two_phases_p
#

def call_stress_strain_two_phases_p(step): # builds stress-strain curve using threads
    pressurefiles=["pressure0", "pressure1", "pressure2", "pressure3"]

    pressure = []

    p0 = threading.Thread(target = stress_strain_two_phases_p, 
                          args = [0, pressurefiles[0], step])
    p1 = threading.Thread(target = stress_strain_two_phases_p,
                          args=[1, pressurefiles[1], step])
    p2 = threading.Thread(target = stress_strain_two_phases_p,
                          args=[2, pressurefiles[2], step])
    p3 = threading.Thread(target = stress_strain_two_phases_p,
                          args=[3, pressurefiles[3], step])

    p0.start()
    p1.start()
    p2.start()
    p3.start()

    for pressurefile in pressurefiles:
        f = open(pressurefile, 'r')
        for line in f:
            pressure.append(line.split())

    for i in range(len(pressure)):
        print(pressure[i][0], pressure[i][1],
              pressure[i][2], pressure[i][3])

    return None

#################################
#
# main functions are called here: 
#

#atoms_displacement() 	# to calculate diffusion coeff curve

#bond_orientation_parameter() # to calculate bond or bond 1-3 
                              # orientation parameter

#phase_distribution() # to get phase distribution

surfactant_distribution() # to estimate surfactant distribution

#stress_strain_two_phases() # builds stress-strain curve

#call_stress_strain_two_phases_p(10) # builds stress-strain curve using threads
