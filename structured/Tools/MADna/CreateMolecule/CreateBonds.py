#!/usr/bin/env python3

import sys
import math
import numpy as np
import glob
import scipy.optimize as optimization
import os
from itertools import islice
import re

#Bond Types
# 1-16 -> SP
# 17-32 -> PS
# 33-36 -> SB
# 37-38 -> BB-inter
# 39-54 -> BB-intra

WC = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

typeBond = {
    'SP/AA': 1,
    'SP/AC': 2,
    'SP/AG': 3,
    'SP/AT': 4,
    'SP/CA': 5,
    'SP/CC': 6,
    'SP/CG': 7,
    'SP/CT': 8,
    'SP/GA': 9,
    'SP/GC': 10,
    'SP/GG': 11,
    'SP/GT': 12,
    'SP/TA': 13,
    'SP/TC': 14,
    'SP/TG': 15,
    'SP/TT': 16,
    'PS/AA': 17,
    'PS/AC': 18,
    'PS/AG': 19,
    'PS/AT': 20,
    'PS/CA': 21,
    'PS/CC': 22,
    'PS/CG': 23,
    'PS/CT': 24,
    'PS/GA': 25,
    'PS/GC': 26,
    'PS/GG': 27,
    'PS/GT': 28,
    'PS/TA': 29,
    'PS/TC': 30,
    'PS/TG': 31,
    'PS/TT': 32,
    'SB/A': 33,
    'SB/C': 34,
    'SB/G': 35,
    'SB/T': 36,
    'BB-inter/AT': 37,
    'BB-inter/CG': 38,
    'BB-inter/GC': 38,
    'BB-inter/TA': 37,
    'BB-intra/AA': 39,
    'BB-intra/AC': 40,
    'BB-intra/AG': 41,
    'BB-intra/AT': 42,
    'BB-intra/CA': 43,
    'BB-intra/CC': 44,
    'BB-intra/CG': 45,
    'BB-intra/CT': 46,
    'BB-intra/GA': 47,
    'BB-intra/GC': 48,
    'BB-intra/GG': 49,
    'BB-intra/GT': 50,
    'BB-intra/TA': 51,
    'BB-intra/TC': 52,
    'BB-intra/TG': 53,
    'BB-intra/TT': 54
}

kbond = {
    1: 108.849600,
    2: 116.093000,
    3: 130.551200,
    4: 134.298200,
    5: 117.795600,
    6: 158.718400,
    7: 121.876200,
    8: 150.068000,
    9: 98.948800,
    10: 82.660600,
    11: 150.242600,
    12: 109.797800,
    13: 139.538000,
    14: 147.736200,
    15: 115.540000,
    16: 159.897200,
    17: 29.178200,
    18: 30.063800,
    19: 23.022000,
    20: 31.743400,
    21: 29.641600,
    22: 38.022200,
    23: 20.852400,
    24: 37.382800,
    25: 30.260400,
    26: 38.421000,
    27: 33.503600,
    28: 34.214400,
    29: 38.532200,
    30: 39.623000,
    31: 32.780400,
    32: 47.345800,
    33: 78.459400,
    34: 89.029200,
    35: 98.264400,
    36: 90.502800,
    37: 42.968000,
    38: 71.908800,
    39: 11.723300,
    40: 16.033040,
    41: 6.953480,
    42: 37.138400,
    43: 11.678180,
    44: 4.556460,
    45: 8.929340,
    46: 8.008620,
    47: 12.141360,
    48: 19.043940,
    49: 6.154040,
    50: 26.376000,
    51: 16.652460,
    52: 5.911920,
    53: 7.927460,
    54: 9.313280
}

r0bond = {
    1: 3.737680,
    2: 3.739850,
    3: 3.750370,
    4: 3.747050,
    5: 3.758460,
    6: 3.765150,
    7: 3.760470,
    8: 3.751930,
    9: 3.738060,
    10: 3.707210,
    11: 3.760770,
    12: 3.741210,
    13: 3.758390,
    14: 3.763810,
    15: 3.758300,
    16: 3.755310,
    17: 4.083580,
    18: 4.070020,
    19: 4.094660,
    20: 4.115280,
    21: 4.127790,
    22: 4.169700,
    23: 4.087590,
    24: 4.119970,
    25: 4.108820,
    26: 4.038130,
    27: 4.173860,
    28: 4.113180,
    29: 4.143020,
    30: 4.166030,
    31: 4.133930,
    32: 4.138290,
    33: 4.894800,
    34: 4.393840,
    35: 5.012620,
    36: 4.457900,
    37: 6.089870,
    38: 5.700030,
    39: 3.871730,
    40: 3.733100,
    41: 4.123320,
    42: 3.696230,
    43: 4.235810,
    44: 4.185290,
    45: 4.259420,
    46: 3.946020,
    47: 3.819480,
    48: 3.708220,
    49: 4.148530,
    50: 3.692320,
    51: 4.358370,
    52: 4.188390,
    53: 4.483310,
    54: 4.016420
}

def create_bonds(seqstring):
	lseq = len(seqstring)
	nbonds = 9*lseq-6
	bondtype = np.zeros(nbonds+1, dtype=int)
	listkbond = np.zeros(nbonds+1, dtype=float)
	listr0bond = np.zeros(nbonds+1, dtype=float)
	i1list = np.zeros(nbonds+1, dtype=int)
	i2list = np.zeros(nbonds+1, dtype=int)
	
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
	
	ibond = 1
	
	#### SP
	for i in range(1,lseq):
	    i1list[ibond] = 3*i - 2
	    i2list[ibond] = i1list[ibond] + 2
	    bondtype[ibond] = typeBond['SP/'+seqstring[i-1]+seqstring[i]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	for i in range(1,lseq):
	    i1list[ibond] = lastB1 + 3*i - 2
	    i2list[ibond] = i1list[ibond] + 2
	    bondtype[ibond] = typeBond['SP/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	
	#### PS
	for i in range(2,lseq+1):
	    i1list[ibond] = 3*i - 3
	    i2list[ibond] = i1list[ibond] + 1
	    bondtype[ibond] = typeBond['PS/'+seqstring[i-2]+seqstring[i-1]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	for i in range(2,lseq+1):
	    i1list[ibond] = lastB1 + 3*i - 3
	    i2list[ibond] = i1list[ibond] + 1
	    bondtype[ibond] = typeBond['PS/'+WC[seqstring[lseq-1-(i-2)]]+WC[seqstring[lseq-1-(i-1)]]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	
	#### SB
	for i in range(1,lseq+1):
	    i1list[ibond] = 3*i - 2
	    i2list[ibond] = i1list[ibond] + 1
	    bondtype[ibond] = typeBond['SB/'+seqstring[i-1]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	for i in range(1,lseq+1):
	    i1list[ibond] = lastB1 + 3*i - 2
	    i2list[ibond] = i1list[ibond] + 1
	    bondtype[ibond] = typeBond['SB/'+WC[seqstring[lseq-1-(i-1)]]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	
	#### BB-inter
	for i in range(1,lseq+1):
	    i1list[ibond] = 3*i - 1
	    i2list[ibond] = lastB2 + 2 - i1list[ibond]
	    bondtype[ibond] = typeBond['BB-inter/'+seqstring[i-1]+WC[seqstring[i-1]]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	
	#### BB-intra
	for i in range(1,lseq):
	    i1list[ibond] = 3*i - 1
	    i2list[ibond] = i1list[ibond] + 3
	    bondtype[ibond] = typeBond['BB-intra/'+seqstring[i-1]+seqstring[i]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	for i in range(1,lseq):
	    i1list[ibond] = lastB1 + 3*i - 1
	    i2list[ibond] = i1list[ibond] + 3
	    bondtype[ibond] = typeBond['BB-intra/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
	    listkbond[ibond] = kbond[bondtype[ibond]]
	    listr0bond[ibond] = r0bond[bondtype[ibond]]
	    ibond += 1
	
	return i1list, i2list, listr0bond, listkbond
