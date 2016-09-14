#!/usr/bin/env python3

import math

######################################################
#
# reads bonds from file only with bonds, returns list:
# [bond_center_z, bond_x, bond_y, bond_z]
#

def read_bonds(datafile):
    f = open(datafile, 'r')
    bonds = [[]]

    # appending bond_center_z and bond vector coordinates
    for string in f:
        string_splitted = string.split()
        bonds.append([int(string_splitted[1]), int(string_splitted[2]),
                      int(string_splitted[3]), 0, 0, 0, 0])     
 
    f.close

    return bonds


##############################################################
#
# returns list (order parameter, number_of bonds) for z layers
#

def bond_orientation(atoms, bonds, step, orientation, mode, scale):
    bounds_mix=[2.3278432446331329e-01,  9.3900097250691388e+01,
               -1.6605997607843719e+00,  7.8477007265866263e+01,
                4.7355184082675912e+00, 4.6138359869416014e+01]
    bounds_seg=[7.5688834338448907e-01, 9.4412591331513482e+01,
               -9.1913280843088074e-01, 7.9300611855381732e+01,
                5.1794498852073787e+01, 9.3240211033719930e+01]

    bounds=bounds_seg

    for i in range(1, len(bonds)):
        bonds[i][3] = (atoms[int(bonds[i][1])][3] +
                       atoms[int(bonds[i][2])][3]) / 2
        dx = abs(atoms[int(bonds[i][1])][1] - atoms[int(bonds[i][2])][1])
        bonds[i][4] = min(dx, abs(bounds[1] - bounds[0] - dx))
        dy = abs(atoms[int(bonds[i][1])][2] - atoms[int(bonds[i][2])][2])
        bonds[i][5] = min(dy, abs(bounds[3] - bounds[2] - dy))
        dz = abs(atoms[int(bonds[i][1])][3] - atoms[int(bonds[i][2])][3])
        bonds[i][6] = min(dz, abs(bounds[5] - bounds[4] - dz))

    if step==1:
        orientation = [[]]
        for i in range(1, len(bonds) * scale):
            orientation.append([0, 0])#number of bonds and order parameter

    for i in range(1, len(bonds)):
        #only main chain bonds
        if mode==1 and (bonds[i][0]!=5 and bonds[i][0]!=8 and
                        bonds[i][0]!=11 and bonds[i][0]!=12 and
                        bonds[i][0]!=13 and bonds[i][0]!=14):
            continue

        #only this phase
        if mode==2 and define_phase(i)!=2 or (bonds[i][0]!=5 and
                                              bonds[i][0]!=8 and 
                                              bonds[i][0]!=11 and
                                              bonds[i][0]!=12 and
                                              bonds[i][0]!=13 and
                                              bonds[i][0]!=14):
            continue
      
        bond_len=math.sqrt(bonds[i][4]**2+bonds[i][5]**2+bonds[i][6]**2)
        cos_alpha=float(bonds[i][6])/bond_len
        if abs(cos_alpha)>1:
            print("ERROR!!!")
        order_par=(3*cos_alpha**2-1)/2
        orientation[int(bonds[i][3]*scale)][0]+=order_par
        orientation[int(bonds[i][3]*scale)][1]+=1

    return orientation
