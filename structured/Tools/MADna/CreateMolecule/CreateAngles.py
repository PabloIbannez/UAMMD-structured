#!/usr/bin/env python3

import sys
import math
import numpy as np
import glob
import scipy.optimize as optimization
import os
from itertools import islice
import re

#Angle Types:
# 1-16 -> SPS
# 17-20 -> SBB
# 21-36 -> 3PSB5
# 37-52 -> 5PSB3

WC = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

typeAngle = {
    'SPS/AA': 1,
    'SPS/AC': 2,
    'SPS/AG': 3,
    'SPS/AT': 4,
    'SPS/CA': 5,
    'SPS/CC': 6,
    'SPS/CG': 7,
    'SPS/CT': 8,
    'SPS/GA': 9,
    'SPS/GC': 10,
    'SPS/GG': 11,
    'SPS/GT': 12,
    'SPS/TA': 13,
    'SPS/TC': 14,
    'SPS/TG': 15,
    'SPS/TT': 16,
    'SBB/AT': 17,
    'SBB/CG': 18,
    'SBB/GC': 19,
    'SBB/TA': 20,
    '3PSB5/AA': 21,
    '3PSB5/AC': 22,
    '3PSB5/AG': 23,
    '3PSB5/AT': 24,
    '3PSB5/CA': 25,
    '3PSB5/CC': 26,
    '3PSB5/CG': 27,
    '3PSB5/CT': 28,
    '3PSB5/GA': 29,
    '3PSB5/GC': 30,
    '3PSB5/GG': 31,
    '3PSB5/GT': 32,
    '3PSB5/TA': 33,
    '3PSB5/TC': 34,
    '3PSB5/TG': 35,
    '3PSB5/TT': 36,
    '5PSB3/AA': 37,
    '5PSB3/AC': 38,
    '5PSB3/AG': 39,
    '5PSB3/AT': 40,
    '5PSB3/CA': 41,
    '5PSB3/CC': 42,
    '5PSB3/CG': 43,
    '5PSB3/CT': 44,
    '5PSB3/GA': 45,
    '5PSB3/GC': 46,
    '5PSB3/GG': 47,
    '5PSB3/GT': 48,
    '5PSB3/TA': 49,
    '5PSB3/TC': 50,
    '5PSB3/TG': 51,
    '5PSB3/TT': 52
}

kangle = {
    1: 52.297400,
    2: 55.454400,
    3: 47.686800,
    4: 66.077800,
    5: 45.379000,
    6: 64.938800,
    7: 43.603800,
    8: 76.369000,
    9: 45.815600,
    10: 50.919600,
    11: 68.554600,
    12: 57.550000,
    13: 52.805600,
    14: 68.361800,
    15: 59.717400,
    16: 93.246000,
    17: 103.155400,
    18: 107.712400,
    19: 104.304600,
    20: 99.360400,
    21: 64.937000,
    22: 81.014400,
    23: 71.817000,
    24: 87.157400,
    25: 76.118000,
    26: 102.752800,
    27: 75.546800,
    28: 102.550000,
    29: 58.532800,
    30: 66.371800,
    31: 105.403600,
    32: 74.350600,
    33: 81.635800,
    34: 97.176400,
    35: 85.435000,
    36: 115.830400,
    37: 12.424360,
    38: 41.473600,
    39: 4.879880,
    40: 40.484000,
    41: 26.925600,
    42: 38.385200,
    43: 16.429200,
    44: 8.518380,
    45: 3.122240,
    46: 20.988400,
    47: 22.892600,
    48: 52.965400,
    49: 17.571740,
    50: 6.500380,
    51: 30.029600,
    52: 26.309600
}

theta0angle = {
    1: 1.642672,
    2: 1.602450,
    3: 1.642474,
    4: 1.620228,
    5: 1.679752,
    6: 1.649765,
    7: 1.666130,
    8: 1.620831,
    9: 1.660466,
    10: 1.617227,
    11: 1.645032,
    12: 1.648561,
    13: 1.647350,
    14: 1.637375,
    15: 1.625242,
    16: 1.617213,
    17: 2.693183,
    18: 2.414680,
    19: 2.770780,
    20: 2.317780,
    21: 2.012836,
    22: 2.026240,
    23: 2.011928,
    24: 2.002451,
    25: 2.080747,
    26: 2.051792,
    27: 2.103087,
    28: 2.045683,
    29: 1.950824,
    30: 1.953687,
    31: 1.913579,
    32: 1.942115,
    33: 2.087745,
    34: 2.112005,
    35: 2.095948,
    36: 2.092440,
    37: 1.975940,
    38: 1.907418,
    39: 2.081096,
    40: 1.846366,
    41: 1.928449,
    42: 1.943127,
    43: 1.978174,
    44: 1.862109,
    45: 1.890907,
    46: 1.861935,
    47: 2.128202,
    48: 1.830554,
    49: 1.994266,
    50: 1.875985,
    51: 2.070519,
    52: 1.831810
}

def create_angles(seqstring):
   lseq = len(seqstring)
   nangles = 8*lseq-6
   angletype = np.zeros(nangles+1, dtype=int)
   i1list = np.zeros(nangles+1, dtype=int)
   i2list = np.zeros(nangles+1, dtype=int)
   i3list = np.zeros(nangles+1, dtype=int)
   listkangle = np.zeros(nangles+1, dtype=float)
   listtheta0angle = np.zeros(nangles+1, dtype=float)
   
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
   
   iangle = 1
   
   #### SPS
   for i in range(1,lseq):
       i1 = 3*i - 2
       i1list[iangle] = i1
       i2list[iangle] = i1 + 2
       i3list[iangle] = i1 + 3
       angletype[iangle] = typeAngle['SPS/'+seqstring[i-1]+seqstring[i]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1
   for i in range(1,lseq):
       i1 = lastB1 + 3*i - 2
       i1list[iangle] = i1
       i2list[iangle] = i1 + 2
       i3list[iangle] = i1 + 3
       angletype[iangle] = typeAngle['SPS/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1
   
   #### SBB
   for i in range(1,lseq+1):
       i1 = 3*i - 2
       i1list[iangle] = i1
       i2list[iangle] = i1 + 1
       i3list[iangle] = lastB2 + 2 - i2list[iangle]
       angletype[iangle] = typeAngle['SBB/'+seqstring[i-1]+WC[seqstring[i-1]]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1
   for i in range(1,lseq+1):
       i1 = lastB1 + 3*i - 2
       i1list[iangle] = i1
       i2list[iangle] = i1 + 1
       i3list[iangle] = lastB2 + 2 - i2list[iangle]
       angletype[iangle] = typeAngle['SBB/'+WC[seqstring[lseq-1-(i-1)]]+seqstring[lseq-1-(i-1)]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1
   
   #### 3PSB5
   for i in range(1,lseq):
       i1 = 3*i
       i1list[iangle] = i1
       i2list[iangle] = i1 - 2
       i3list[iangle] = i1 - 1
       angletype[iangle] = typeAngle['3PSB5/'+seqstring[i-1]+seqstring[i]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i1list[iangle] = i1
       i2list[iangle] = i1 - 2
       i3list[iangle] = i1 - 1
       angletype[iangle] = typeAngle['3PSB5/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1
   
   #### 5PSB3
   for i in range(1,lseq):
       i1 = 3*i
       i1list[iangle] = i1
       i2list[iangle] = i1 + 1
       i3list[iangle] = i1 + 2
       angletype[iangle] = typeAngle['5PSB3/'+seqstring[i-1]+seqstring[i]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i1list[iangle] = i1
       i2list[iangle] = i1 + 1
       i3list[iangle] = i1 + 2
       angletype[iangle] = typeAngle['5PSB3/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       listkangle[iangle] = kangle[angletype[iangle]]
       listtheta0angle[iangle] = theta0angle[angletype[iangle]]
       iangle += 1

   return i1list, i2list, i3list, listtheta0angle, listkangle
