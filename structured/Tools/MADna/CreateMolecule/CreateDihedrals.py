#!/usr/bin/env python3

import sys
import math
import numpy as np
import glob
import scipy.optimize as optimization
import os
from itertools import islice
import re

#Dihedral Types:
# 1-16 -> SPSP
# 17-32 -> PSPS
# 33-48 -> SPSB53
# 49-64 -> SPSB35
# 65-80 -> PSBB53
# 81-96 -> PSBB35

WC = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

typeDihedral = {
    'SPSP/AA': 1,
    'SPSP/AC': 2,
    'SPSP/AG': 3,
    'SPSP/AT': 4,
    'SPSP/CA': 5,
    'SPSP/CC': 6,
    'SPSP/CG': 7,
    'SPSP/CT': 8,
    'SPSP/GA': 9,
    'SPSP/GC': 10,
    'SPSP/GG': 11,
    'SPSP/GT': 12,
    'SPSP/TA': 13,
    'SPSP/TC': 14,
    'SPSP/TG': 15,
    'SPSP/TT': 16,
    'PSPS/AA': 17,
    'PSPS/AC': 18,
    'PSPS/AG': 19,
    'PSPS/AT': 20,
    'PSPS/CA': 21,
    'PSPS/CC': 22,
    'PSPS/CG': 23,
    'PSPS/CT': 24,
    'PSPS/GA': 25,
    'PSPS/GC': 26,
    'PSPS/GG': 27,
    'PSPS/GT': 28,
    'PSPS/TA': 29,
    'PSPS/TC': 30,
    'PSPS/TG': 31,
    'PSPS/TT': 32,
    'SPSB53/AA': 33,
    'SPSB53/AC': 34,
    'SPSB53/AG': 35,
    'SPSB53/AT': 36,
    'SPSB53/CA': 37,
    'SPSB53/CC': 38,
    'SPSB53/CG': 39,
    'SPSB53/CT': 40,
    'SPSB53/GA': 41,
    'SPSB53/GC': 42,
    'SPSB53/GG': 43,
    'SPSB53/GT': 44,
    'SPSB53/TA': 45,
    'SPSB53/TC': 46,
    'SPSB53/TG': 47,
    'SPSB53/TT': 48,
    'SPSB35/AA': 49,
    'SPSB35/AC': 50,
    'SPSB35/AG': 51,
    'SPSB35/AT': 52,
    'SPSB35/CA': 53,
    'SPSB35/CC': 54,
    'SPSB35/CG': 55,
    'SPSB35/CT': 56,
    'SPSB35/GA': 57,
    'SPSB35/GC': 58,
    'SPSB35/GG': 59,
    'SPSB35/GT': 60,
    'SPSB35/TA': 61,
    'SPSB35/TC': 62,
    'SPSB35/TG': 63,
    'SPSB35/TT': 64,
    'PSBB53/AA': 65,
    'PSBB53/AC': 66,
    'PSBB53/AG': 67,
    'PSBB53/AT': 68,
    'PSBB53/CA': 69,
    'PSBB53/CC': 70,
    'PSBB53/CG': 71,
    'PSBB53/CT': 72,
    'PSBB53/GA': 73,
    'PSBB53/GC': 74,
    'PSBB53/GG': 75,
    'PSBB53/GT': 76,
    'PSBB53/TA': 77,
    'PSBB53/TC': 78,
    'PSBB53/TG': 79,
    'PSBB53/TT': 80,
    'PSBB35/AA': 81,
    'PSBB35/AC': 82,
    'PSBB35/AG': 83,
    'PSBB35/AT': 84,
    'PSBB35/CA': 85,
    'PSBB35/CC': 86,
    'PSBB35/CG': 87,
    'PSBB35/CT': 88,
    'PSBB35/GA': 89,
    'PSBB35/GC': 90,
    'PSBB35/GG': 91,
    'PSBB35/GT': 92,
    'PSBB35/TA': 93,
    'PSBB35/TC': 94,
    'PSBB35/TG': 95,
    'PSBB35/TT': 96
}

kdihedral = {
    1: 59.550400,
    2: 23.660000,
    3: 47.340200,
    4: 19.102700,
    5: 25.628000,
    6: 44.821700,
    7: 54.620000,
    8: 57.425500,
    9: 50.067500,
    10: 52.403800,
    11: 42.783600,
    12: 25.742000,
    13: 19.632700,
    14: 60.413500,
    15: 32.650300,
    16: 77.554100,
    17: 4.724600,
    18: 4.821760,
    19: 5.557370,
    20: 7.502810,
    21: 5.319440,
    22: 11.514500,
    23: 4.110620,
    24: 10.462500,
    25: 4.598920,
    26: 2.748350,
    27: 10.890000,
    28: 5.608680,
    29: 8.347130,
    30: 10.759600,
    31: 7.747140,
    32: 15.657500,
    33: 13.817700,
    34: 19.093900,
    35: 10.119900,
    36: 27.135700,
    37: 10.957000,
    38: 21.610000,
    39: 7.563410,
    40: 27.527200,
    41: 16.593400,
    42: 21.405900,
    43: 13.633500,
    44: 28.047400,
    45: 14.212100,
    46: 23.547500,
    47: 9.022560,
    48: 34.024000,
    49: 10.394900,
    50: 10.911200,
    51: 8.938090,
    52: 14.860400,
    53: 6.912550,
    54: 13.923100,
    55: 5.558970,
    56: 13.430800,
    57: 11.527900,
    58: 11.038800,
    59: 15.044200,
    60: 14.542200,
    61: 8.691430,
    62: 12.006000,
    63: 6.706810,
    64: 15.723900,
    65: 10.495200,
    66: 13.636700,
    67: 6.697630,
    68: 23.898700,
    69: 7.624270,
    70: 17.910700,
    71: 6.389590,
    72: 18.040600,
    73: 10.672300,
    74: 15.305400,
    75: 7.243500,
    76: 24.122700,
    77: 8.404690,
    78: 12.935000,
    79: 6.903510,
    80: 21.298700,
    81: 10.116000,
    82: 7.446700,
    83: 8.558080,
    84: 8.176150,
    85: 8.981520,
    86: 13.919100,
    87: 8.851580,
    88: 10.039300,
    89: 7.921290,
    90: 6.712400,
    91: 7.913630,
    92: 7.556260,
    93: 12.113600,
    94: 12.990100,
    95: 10.703100,
    96: 14.626600
}

phi0dihedral = {
    1: 0.007487,
    2: -0.026250,
    3: -0.037751,
    4: -0.015237,
    5: -0.046094,
    6: -0.064420,
    7: -0.013840,
    8: -0.008779,
    9: -0.008290,
    10: 0.017087,
    11: -0.055275,
    12: 0.012235,
    13: -0.049882,
    14: 0.001187,
    15: -0.046513,
    16: 0.035256,
    17: 0.390256,
    18: 0.478081,
    19: 0.466509,
    20: 0.426262,
    21: 0.349118,
    22: 0.412037,
    23: 0.463333,
    24: 0.430782,
    25: 0.334998,
    26: 0.309045,
    27: 0.410938,
    28: 0.332189,
    29: 0.385840,
    30: 0.337145,
    31: 0.462582,
    32: 0.358107,
    33: -2.484267,
    34: -2.465906,
    35: -2.544498,
    36: -2.398466,
    37: -2.524340,
    38: -2.449762,
    39: -2.610280,
    40: -2.386179,
    41: -2.440145,
    42: -2.407525,
    43: -2.461211,
    44: -2.359650,
    45: -2.510761,
    46: -2.408816,
    47: -2.564290,
    48: -2.366841,
    49: 2.861886,
    50: 2.901139,
    51: 2.881416,
    52: 2.870229,
    53: 2.716099,
    54: 2.746590,
    55: 2.786087,
    56: 2.774218,
    57: 2.886722,
    58: 2.893686,
    59: 2.909987,
    60: 2.881957,
    61: 2.655117,
    62: 2.606911,
    63: 2.698646,
    64: 2.650754,
    65: 0.365402,
    66: 0.345436,
    67: 0.391460,
    68: 0.374146,
    69: 0.403642,
    70: 0.150325,
    71: 0.667152,
    72: 0.373919,
    73: 0.370115,
    74: 0.416488,
    75: 0.180607,
    76: 0.388859,
    77: 0.484451,
    78: 0.345994,
    79: 0.626381,
    80: 0.373989,
    81: -2.177857,
    82: -2.098095,
    83: -2.097292,
    84: -2.027113,
    85: -2.042192,
    86: -2.205799,
    87: -1.973706,
    88: -2.093051,
    89: -2.177944,
    90: -1.991351,
    91: -2.278999,
    92: -1.980530,
    93: -1.968068,
    94: -1.957422,
    95: -1.931504,
    96: -2.021266
}

def create_dihedrals(seqstring):
   lseq = len(seqstring)
   ndihedrals = 12*lseq-16
   dihedraltype = np.zeros(ndihedrals+1, dtype=int)
   i1list = np.zeros(ndihedrals+1, dtype=int)
   i2list = np.zeros(ndihedrals+1, dtype=int)
   i3list = np.zeros(ndihedrals+1, dtype=int)
   i4list = np.zeros(ndihedrals+1, dtype=int)
   listkdihedral = np.zeros(ndihedrals+1, dtype=float)
   listphi0dihedral = np.zeros(ndihedrals+1, dtype=float)
   
   firstS1 = 1
   lastS1 = firstS1 + (lseq-1)*3
   firstP1 = 3
   lastP1 = firstP1 + (lseq-2)*3
   firstB1 = 2
   lastB1 = firstB1 + (lseq-1)*3
   
   firstS2 = lastB1 + 1
   lastS2 = firstS2 + (lseq-1)*3
   firstP2 = lastB1 + 3
   lastP2 = firstP2 + (lseq-2)*3
   firstB2 = lastB1 + 2
   lastB2 = firstB2 + (lseq-1)*3
   
   idihedral = 1
   
   #### SPSP
   for i in range(1,lseq-1):
       i1 = 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 5
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['SPSP/'+seqstring[i-1]+seqstring[i]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   for i in range(1,lseq-1):
       i1 = lastB1 + 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 5
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['SPSP/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   
   #### PSPS
   for i in range(2,lseq):
       i1 = 3*i - 3
       i2 = i1 + 1
       i3 = i1 + 3
       i4 = i1 + 4
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['PSPS/'+seqstring[i-1]+seqstring[i]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   for i in range(2,lseq):
       i1 = lastB1 + 3*i - 3
       i2 = i1 + 1
       i3 = i1 + 3
       i4 = i1 + 4
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['PSPS/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   
   #### SPSB53
   for i in range(1,lseq):
       i1 = 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 4
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['SPSB53/'+seqstring[i-1]+seqstring[i]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   for i in range(1,lseq):
       i1 = lastB1 + 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 4
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['SPSB53/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   
   #### SPSB35
   for i in range(1,lseq):
       i1 = 3*i - 1
       i2 = i1 - 1
       i3 = i1 + 1
       i4 = i1 + 2
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['SPSB35/'+seqstring[i-1]+seqstring[i]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   for i in range(1,lseq):
       i1 = lastB1 + 3*i - 1
       i2 = i1 - 1
       i3 = i1 + 1
       i4 = i1 + 2
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['SPSB35/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   
   #### PSBB53
   for i in range(1,lseq):
       i1 = 3*i
       i2 = i1 + 1
       i3 = i1 + 2
       i4 = lastB2 + 2 - i3
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['PSBB53/'+seqstring[i-1]+seqstring[i]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i2 = i1 + 1
       i3 = i1 + 2
       i4 = lastB2 + 2 - i3
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['PSBB53/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   
   #### PSBB35
   for i in range(1,lseq):
       i1 = 3*i
       i2 = i1 - 2
       i3 = i1 - 1
       i4 = lastB2 + 2 - i3
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['PSBB35/'+seqstring[i-1]+seqstring[i]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i2 = i1 - 2
       i3 = i1 - 1
       i4 = lastB2 + 2 - i3
       i1list[idihedral] = i1
       i2list[idihedral] = i2
       i3list[idihedral] = i3
       i4list[idihedral] = i4
       dihedraltype[idihedral] = typeDihedral['PSBB35/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkdihedral[idihedral] = kdihedral[dihedraltype[idihedral]]
       listphi0dihedral[idihedral] = phi0dihedral[dihedraltype[idihedral]]
       idihedral += 1
   
   return i1list, i2list, i3list, i4list, listphi0dihedral, listkdihedral
