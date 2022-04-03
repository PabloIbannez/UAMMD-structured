#!/usr/bin/env python3

import sys
import os
import math
import numpy as np
import copy
import glob
from itertools import islice
import re

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

scr_fold = get_script_path()
sys.path.insert(1, scr_fold)

import CreateCoordinates
import CreateBonds
import CreateAngles
import CreateDihedrals

if (len(sys.argv) != 3+1):
    print("Usage:") 
    print("# 1 -> input sequence")
    print("# 2 -> initial coordinate file")
    print("# 3 -> topology file")
    sys.exit(0)

seqstring = sys.argv[1]
datafile = sys.argv[2]
topfile = sys.argv[3]

file_writer=open(topfile,'w')

def stampa(*argv):
    s = ""
    for arg in argv:
        s += str(arg)+" "
    print(s, file=file_writer)

### parameters
lseq = len(seqstring)
natoms = lseq*6-2
nbonds = lseq*9-6
nangles = lseq*8-6
ndihedrals = lseq*12-16

### create coordinates and topology
coordlist, typelist, chargelist = CreateCoordinates.create_coordinates(seqstring, datafile)
i1bond, i2bond, r0bond, kbond = CreateBonds.create_bonds(seqstring)
i1angle, i2angle, i3angle, theta0angle, kangle = CreateAngles.create_angles(seqstring)
i1dihedral, i2dihedral, i3dihedral, i4dihedral, phi0dihedral, kdihedral = CreateDihedrals.create_dihedrals(seqstring)

### dictionary connecting type and name
tipo = {
    1: 'S',
    2: 'P',
    3: 'A',
    4: 'C',
    5: 'G',
    6: 'T'
}

### create exclusions matrix
exclusions = []
for i in range(0, natoms+1):
        exclusions.append([])
neighbors = []
for i in range(0, natoms+1):
        neighbors.append([])
for ibond in range(1, nbonds+1):
        i1 = i1bond[ibond]
        i2 = i2bond[ibond]
        neighbors[i1].append(i2)
        neighbors[i2].append(i1)
for i in range(1, natoms+1):
        vicini_actual = [i]
        for distanza in range(1,4):
                vicini_next = []
                for i0 in vicini_actual:
                        for j in neighbors[i0]:
                                vicini_next.append(j)
                                exclusions[i].append(j)
                vicini_actual = vicini_next
for i in range(0, natoms+1):
        exclusions[i] = list(dict.fromkeys(exclusions[i]))
        exclusions[i].sort(key=int)

### print topology file
stampa("[TYPES]")
stampa("S 76.05 3.2 0.0")
stampa("P 94.97 2.25 -0.6")
stampa("A 130.09 2.7 0.0")
stampa("C 106.06 3.2 0.0")
stampa("G 146.09 2.45 0.0")
stampa("T 120.07 3.55 0.0")

stampa("[STRUCTURE]")
nucID = 0
for i in range(1,natoms+1):
    if i <= natoms/2:
        chainID = 0
    else:
        chainID = 1
    stampa(i-1, tipo[typelist[i]], nucID, chainID, 0)
    if tipo[typelist[i]] != "S" and tipo[typelist[i]] != "P":
        nucID = nucID + 1

stampa("[BONDS]")
for i in range(1, nbonds+1):
	stampa(i1bond[i]-1, i2bond[i]-1, r0bond[i], kbond[i])

stampa("[ANGLES]")
for i in range(1, nangles+1):
	stampa(i1angle[i]-1, i2angle[i]-1, i3angle[i]-1, theta0angle[i], kangle[i])

stampa("[DIHEDRALS]")
for i in range(1, ndihedrals+1):
	stampa(i1dihedral[i]-1, i2dihedral[i]-1, i3dihedral[i]-1, i4dihedral[i]-1, 1, phi0dihedral[i], kdihedral[i])

stampa("[EXCLUSIONS]")
for i in range(1, natoms+1):
    stringa = str(i-1) + " " 
    for j in exclusions[i]:
        if j != i:
            stringa += str(j-1) + " "
    stampa(stringa)
