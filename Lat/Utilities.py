#!/usr/bin/env python3

################################################
#
# depending on atom number returns phase number:
# 1 for mmt
# 2 for modifier
# 3 for polymer
#

def define_phase(i):

#----------------mmt
    if (i < 721
    or (i > 3480 and i < 4201)
    or (i > 6960 and i < 7681)
    or (i > 10440 and i < 11160)
    or (i > 13920 and i < 14641)
    or (i > 17420 and i < 18121)
    or (i > 20880 and i < 21601)
    or (i > 24360 and i < 25081)
    or (i > 27840 and i < 28561)):
        return 1 
   
#---------------modifier
    elif ((i > 720 and i < 1561)
    or (i > 4200 and i < 5041)
    or (i > 7680 and i < 8521)
    or (i > 11160 and i < 12001)
    or (i > 14640 and i < 15481)
    or (i > 18120 and i < 18961)
    or (i > 21600 and i < 22441)
    or (i > 25080 and i < 25921)
    or (i > 28560 and i < 29401)):
        return 2

#---------------polymer
    else:
        return 3


#######################################################################
#
# calculates center of mass for atoms array, returns [cm_x, cm_y, cm_z]
#

def calculate_cm(atoms):
    cm=[0, 0, 0]	#cm coordinates
    mass=0		#whole mass

    masses=[26.9815, 24.305, 28.0855, 15.9994, 15.9994, 15.9994, 15.9994,
           1.00797, 14.0067, 12.0112, 12.0112, 1.00797, 1.00797, 12.0112,
           15.9994, 14.0067, 14.0067]  #masses from lammps datafile

    for i in range (1, len(atoms)):
        cm[0] += atoms[i][1] * masses[int(atoms[i][0]) - 1]
        cm[1] += atoms[i][2] * masses[int(atoms[i][0]) - 1]
        cm[2] += atoms[i][3] * masses[int(atoms[i][0]) - 1]
        mass += masses[int(atoms[i][0]) - 1]

    cm[0] /= mass
    cm[1] /= mass
    cm[2] /= mass

    return cm
