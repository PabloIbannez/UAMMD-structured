#!/usr/bin/env python3

import sys
import math
import numpy as np
import copy
import glob
import os
from itertools import islice
import re

def stampa(*argv):
    s = ""
    for arg in argv:
        s += str(arg)+" "
    print(s)

def acos(x):
    if x > 1:
        x = 1
    if x < -1:
        x = -1
    return np.arccos(x)

def align(ref_snap, traj_snap):
    centro_ref = np.mean(ref_snap, axis=0)
    centro_traj = np.mean(traj_snap, axis=0)

    ######## translate coordinates so that centroid is 0 
    ref_snap = ref_snap - centro_ref
    traj_snap = traj_snap - centro_traj

    ######## Kabsch algorithm
    ## 1) compute covariance matrix C
    C = np.dot(np.transpose(traj_snap), ref_snap)
    ## 2) compute singular value decomposition
    V, S, W = np.linalg.svd(C)
    ## 3) eventually correct sign to ensure right-handedness
    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    ## 4) compute rotation matrix U
    U = np.dot(V, W)

    ######## rotate traj_snap and translate it onto the original centroid of ref
    traj_snap = np.dot(traj_snap, U)
    traj_snap = traj_snap + centro_ref
    ref_snap = ref_snap + centro_ref

    ######## compute rmsd
    diff = traj_snap - ref_snap
    rmsd = np.sqrt((diff*diff).sum()/len(traj_snap))

    return rmsd, U, traj_snap

WC = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

tipo = {
    'S': 1,
    'P': 2,
    'A': 3,
    'C': 4,
    'G': 5,
    'T': 6
}

charge = {
    'S': 0.0,
    'P': -0.6,
    'A': 0.0,
    'C': 0.0,
    'G': 0.0,
    'T': 0.0
    }

### S,B,P stand for Sugar, Base, Phosphate
### 5,3 refers to the position of the particle within the strand (reference is the 5'3' direction of strand 1)
### 1,2 is the strand

pos_AA = {
    'S_5_1': np.array([104.291942, 11.354422, 44.665789]),
    'B_5_1': np.array([108.186569, 14.007470, 43.378059]),
    'B_5_2': np.array([114.135259, 15.272054, 43.205896]),
    'S_5_2': np.array([117.576789, 13.957587, 45.693125]),
    'P_1': np.array([103.847734, 10.795693, 48.317319]),
    'P_2': np.array([120.511591, 16.787221, 45.117158]),
    'S_3_1': np.array([107.642218, 11.848946, 49.262477]),
    'B_3_1': np.array([109.393357, 15.597971, 46.667270]),
    'B_3_2': np.array([113.829391, 19.074166, 44.373470]),
    'S_3_2': np.array([118.103149, 19.626170, 45.469636])
}

pos_AC = {
    'S_5_1': np.array([84.053976, -102.727702, -12.419664]),
    'B_5_1': np.array([80.271694, -104.425330, -9.835851]),
    'B_5_2': np.array([77.458552, -108.386836, -6.173238]),
    'S_5_2': np.array([78.173519, -112.661152, -5.198133]),
    'P_1': np.array([85.306732, -104.875691, -15.188228]),
    'P_2': np.array([74.574626, -114.052808, -3.823005]),
    'S_3_1': np.array([83.292734, -108.112922, -13.840258]),
    'B_3_1': np.array([79.325730, -106.948985, -12.396229]),
    'B_3_2': np.array([75.251405, -109.051833, -9.021899]),
    'S_3_2': np.array([72.579734, -112.480740, -6.546104])
}

pos_AG = {
    'S_5_1': np.array([104.478177, 157.643352, 207.881915]),
    'B_5_1': np.array([107.752533, 157.793503, 204.261831]),
    'B_5_2': np.array([110.433858, 155.694413, 199.221458]),
    'S_5_2': np.array([110.507018, 151.694831, 197.288646]),
    'P_1': np.array([104.681724, 154.671944, 210.126573]),
    'P_2': np.array([113.757850, 151.558448, 194.808018]),
    'S_3_1': np.array([107.009577, 152.522444, 207.600031]),
    'B_3_1': np.array([110.470242, 154.774466, 204.778691]),
    'B_3_2': np.array([114.130254, 155.545130, 200.486414]),
    'S_3_2': np.array([116.004765, 153.179470, 197.311423])
}

pos_AT = {
    'S_5_1': np.array([48.920256, 72.737169, -153.332677]),
    'B_5_1': np.array([46.610067, 69.414506, -150.596312]),
    'B_5_2': np.array([42.579094, 64.962232, -149.622468]),
    'S_5_2': np.array([38.842468, 63.760005, -151.700334]),
    'P_1': np.array([46.781074, 75.520778, -154.587645]),
    'P_2': np.array([36.676721, 61.558924, -149.009742]),
    'S_3_1': np.array([43.299642, 73.523946, -153.769200]),
    'B_3_1': np.array([43.826045, 71.662836, -149.770242]),
    'B_3_2': np.array([41.076668, 67.012179, -146.972635]),
    'S_3_2': np.array([37.389665, 63.901826, -146.198746])
}

pos_CA = {
    'S_5_1': np.array([83.066905, -108.135394, -14.233104]),
    'B_5_1': np.array([79.310752, -106.952640, -12.312937]),
    'B_5_2': np.array([75.315419, -109.303574, -9.010190]),
    'S_5_2': np.array([72.543516, -112.899742, -6.914360]),
    'P_1': np.array([83.050555, -110.817541, -16.835136]),
    'P_2': np.array([68.612370, -112.170594, -5.993816]),
    'S_3_1': np.array([79.657897, -112.699786, -15.538350]),
    'B_3_1': np.array([76.357327, -109.473991, -13.943321]),
    'B_3_2': np.array([71.516654, -107.710230, -10.708287]),
    'S_3_2': np.array([68.091606, -109.910509, -8.923989])
}

pos_CC = {
    'S_5_1': np.array([159.447965, 54.266892, -40.497196]),
    'B_5_1': np.array([155.704740, 56.509640, -40.086452]),
    'B_5_2': np.array([153.865267, 61.864574, -40.686060]),
    'S_5_2': np.array([153.996776, 66.734806, -41.818849]),
    'P_1': np.array([162.438599, 55.510826, -38.611195]),
    'P_2': np.array([150.959728, 69.376056, -40.800944]),
    'S_3_1': np.array([161.051008, 59.408327, -38.285344]),
    'B_3_1': np.array([156.758306, 59.543681, -37.429257]),
    'B_3_2': np.array([152.600504, 63.419528, -37.098128]),
    'S_3_2': np.array([150.356107, 67.878271, -37.422676])
}

pos_CG = {
    'S_5_1': np.array([107.175975, -66.720492, 116.382896]),
    'B_5_1': np.array([106.529167, -66.466159, 120.708297]),
    'B_5_2': np.array([102.606525, -68.032563, 124.524953]),
    'S_5_2': np.array([98.094046, -69.104042, 126.391554]),
    'P_1': np.array([105.051110, -64.755742, 114.015529]),
    'P_2': np.array([97.520950, -67.973509, 130.233407]),
    'S_3_1': np.array([101.951128, -64.242316, 116.561693]),
    'B_3_1': np.array([103.215332, -63.915418, 121.387619]),
    'B_3_2': np.array([102.845358, -64.617855, 127.024860]),
    'S_3_2': np.array([99.729310, -64.960504, 130.085193])
}

pos_CT = {
    'S_5_1': np.array([116.004765, 153.179470, 197.311423]),
    'B_5_1': np.array([114.130254, 155.545130, 200.486414]),
    'B_5_2': np.array([110.470242, 154.774466, 204.778691]),
    'S_5_2': np.array([107.009577, 152.522444, 207.600031]),
    'P_1': np.array([113.757850, 151.558448, 194.808018]),
    'P_2': np.array([104.681724, 154.671944, 210.126573]),
    'S_3_1': np.array([110.507018, 151.694831, 197.288646]),
    'B_3_1': np.array([110.433858, 155.694413, 199.221458]),
    'B_3_2': np.array([107.752533, 157.793503, 204.261831]),
    'S_3_2': np.array([104.478177, 157.643352, 207.881915])
}

pos_GA = {
    'S_5_1': np.array([106.807165, 152.356498, 207.144812]),
    'B_5_1': np.array([110.444759, 154.631521, 204.574604]),
    'B_5_2': np.array([114.254558, 155.517622, 200.437530]),
    'S_5_2': np.array([116.107197, 153.293915, 197.151307]),
    'P_1': np.array([107.704742, 148.876557, 208.098304]),
    'P_2': np.array([119.574958, 155.270727, 196.079168]),
    'S_3_1': np.array([110.847786, 148.558477, 205.541186]),
    'B_3_1': np.array([113.547739, 152.622611, 205.363927]),
    'B_3_2': np.array([117.541192, 156.434313, 202.809707]),
    'S_3_2': np.array([120.595904, 155.905759, 199.625456])
}

pos_GC = {
    'S_5_1': np.array([102.090213, -64.293508, 116.696675]),
    'B_5_1': np.array([103.362385, -64.035820, 121.526541]),
    'B_5_2': np.array([102.958384, -64.754081, 127.160729]),
    'S_5_2': np.array([99.830731, -65.002848, 130.216666]),
    'P_1': np.array([99.092021, -62.300201, 115.877484]),
    'P_2': np.array([101.463618, -64.011089, 133.741677]),
    'S_3_1': np.array([98.036215, -61.796889, 119.712987]),
    'B_3_1': np.array([101.647670, -60.822798, 121.991164]),
    'B_3_2': np.array([103.508726, -61.130240, 127.362865]),
    'S_3_2': np.array([103.324136, -61.135325, 132.360213])
}

pos_GG = {
    'S_5_1': np.array([150.356107, 67.878271, -37.422676]),
    'B_5_1': np.array([152.600504, 63.419528, -37.098128]),
    'B_5_2': np.array([156.758306, 59.543681, -37.429257]),
    'S_5_2': np.array([161.051008, 59.408327, -38.285344]),
    'P_1': np.array([150.959728, 69.376056, -40.800944]),
    'P_2': np.array([162.438599, 55.510826, -38.611195]),
    'S_3_1': np.array([153.996776, 66.734806, -41.818849]),
    'B_3_1': np.array([153.865267, 61.864574, -40.686060]),
    'B_3_2': np.array([155.704740, 56.509640, -40.086452]),
    'S_3_2': np.array([159.447965, 54.266892, -40.497196])
}

pos_GT = {
    'S_5_1': np.array([72.579734, -112.480740, -6.546104]),
    'B_5_1': np.array([75.251405, -109.051833, -9.021899]),
    'B_5_2': np.array([79.325730, -106.948985, -12.396229]),
    'S_5_2': np.array([83.292734, -108.112922, -13.840258]),
    'P_1': np.array([74.574626, -114.052808, -3.823005]),
    'P_2': np.array([85.306732, -104.875691, -15.188228]),
    'S_3_1': np.array([78.173519, -112.661152, -5.198133]),
    'B_3_1': np.array([77.458552, -108.386836, -6.173238]),
    'B_3_2': np.array([80.271694, -104.425330, -9.835851]),
    'S_3_2': np.array([84.053976, -102.727702, -12.419664])
}

pos_TA = {
    'S_5_1': np.array([43.275101, 73.865246, -153.545926]),
    'B_5_1': np.array([43.645123, 71.853394, -149.598850]),
    'B_5_2': np.array([40.761336, 67.123208, -147.086356]),
    'S_5_2': np.array([37.030824, 64.015686, -146.563577]),
    'P_1': np.array([40.168114, 75.764523, -154.389521]),
    'P_2': np.array([36.653184, 62.366742, -142.820042]),
    'S_3_1': np.array([37.725916, 73.502627, -151.980324]),
    'B_3_1': np.array([39.734424, 72.081418, -147.761773]),
    'B_3_2': np.array([40.390025, 68.392589, -142.969527]),
    'S_3_2': np.array([37.666653, 65.475667, -141.011104])
}

pos_TC = {
    'S_5_1': np.array([120.595904, 155.905759, 199.625456]),
    'B_5_1': np.array([117.541192, 156.434313, 202.809707]),
    'B_5_2': np.array([113.547739, 152.622611, 205.363927]),
    'S_5_2': np.array([110.847786, 148.558477, 205.541186]),
    'P_1': np.array([119.574958, 155.270727, 196.079168]),
    'P_2': np.array([107.704742, 148.876557, 208.098304]),
    'S_3_1': np.array([116.107197, 153.293915, 197.151307]),
    'B_3_1': np.array([114.254558, 155.517622, 200.437530]),
    'B_3_2': np.array([110.444759, 154.631521, 204.574604]),
    'S_3_2': np.array([106.807165, 152.356498, 207.144812])
}

pos_TG = {
    'S_5_1': np.array([68.091606, -109.910509, -8.923989]),
    'B_5_1': np.array([71.516654, -107.710230, -10.708287]),
    'B_5_2': np.array([76.357327, -109.473991, -13.943321]),
    'S_5_2': np.array([79.657897, -112.699786, -15.538350]),
    'P_1': np.array([68.612370, -112.170594, -5.993816]),
    'P_2': np.array([83.050555, -110.817541, -16.835136]),
    'S_3_1': np.array([72.543516, -112.899742, -6.914360]),
    'B_3_1': np.array([75.315419, -109.303574, -9.010190]),
    'B_3_2': np.array([79.310752, -106.952640, -12.312937]),
    'S_3_2': np.array([83.066905, -108.135394, -14.233104])
}

pos_TT = {
    'S_5_1': np.array([118.103149, 19.626170, 45.469636]),
    'B_5_1': np.array([113.829391, 19.074166, 44.373470]),
    'B_5_2': np.array([109.393357, 15.597971, 46.667270]),
    'S_5_2': np.array([107.642218, 11.848946, 49.262477]),
    'P_1': np.array([120.511591, 16.787221, 45.117158]),
    'P_2': np.array([103.847734, 10.795693, 48.317319]),
    'S_3_1': np.array([117.576789, 13.957587, 45.693125]),
    'B_3_1': np.array([114.135259, 15.272054, 43.205896]),
    'B_3_2': np.array([108.186569, 14.007470, 43.378059]),
    'S_3_2': np.array([104.291942, 11.354422, 44.665789])
}

pos_dic = {
    'AA': pos_AA,
    'AC': pos_AC,
    'AG': pos_AG,
    'AT': pos_AT,
    'CA': pos_CA,
    'CC': pos_CC,
    'CG': pos_CG,
    'CT': pos_CT,
    'GA': pos_GA,
    'GC': pos_GC,
    'GG': pos_GG,
    'GT': pos_GT,
    'TA': pos_TA,
    'TC': pos_TC,
    'TG': pos_TG,
    'TT': pos_TT
}

def update_coordinates(loc, coordlist, typelist, chargelist, iseq, seqstring, isfirst):
    natoms = (len(seqstring)-1)*6 - 2
    iS_5_1 = iseq*3 - 2
    iB_5_1 = iS_5_1 + 1
    iB_5_2 = natoms + 1 - iS_5_1
    iS_5_2 = natoms - iS_5_1
    iP_1 = iS_5_1 + 2
    iP_2 = natoms - iS_5_1 - 1
    iS_3_1 = iS_5_1 + 3
    iB_3_1 = iS_5_1 + 4
    iB_3_2 = natoms - iS_5_1 - 2
    iS_3_2 = natoms - iS_5_1 - 3
    if isfirst == 'y':
        coordlist[iS_5_1] = loc["S_5_1"]
        coordlist[iB_5_1] = loc["B_5_1"]
        coordlist[iB_5_2] = loc["B_5_2"]
        coordlist[iS_5_2] = loc["S_5_2"]
        typelist[iS_5_1] = tipo["S"]
        typelist[iB_5_1] = tipo[seqstring[iseq]]
        typelist[iB_5_2] = tipo[WC[seqstring[iseq]]]
        typelist[iS_5_2] = tipo["S"]
        chargelist[iS_5_1] = charge["S"]
        chargelist[iB_5_1] = charge[seqstring[iseq]]
        chargelist[iB_5_2] = charge[WC[seqstring[iseq]]]
        chargelist[iS_5_2] = charge["S"]
    else:
        coordlist[iS_5_1] = 0.5*(loc["S_5_1"] + coordlist[iS_5_1])
        coordlist[iB_5_1] = 0.5*(loc["B_5_1"] + coordlist[iB_5_1]) 
        coordlist[iB_5_2] = 0.5*(loc["B_5_2"] + coordlist[iB_5_2]) 
        coordlist[iS_5_2] = 0.5*(loc["S_5_2"] + coordlist[iS_5_2]) 
        
    coordlist[iP_1] = loc["P_1"]
    coordlist[iP_2] = loc["P_2"]
    coordlist[iS_3_1] = loc["S_3_1"]
    coordlist[iB_3_1] = loc["B_3_1"]
    coordlist[iB_3_2] = loc["B_3_2"]
    coordlist[iS_3_2] = loc["S_3_2"]
    typelist[iP_1] = tipo["P"]
    typelist[iP_2] = tipo["P"]
    typelist[iS_3_1] = tipo["S"]
    typelist[iB_3_1] = tipo[seqstring[iseq+1]]
    typelist[iB_3_2] = tipo[WC[seqstring[iseq+1]]]
    typelist[iS_3_2] = tipo["S"]
    chargelist[iP_1] = charge["P"]
    chargelist[iP_2] = charge["P"]
    chargelist[iS_3_1] = charge["S"]
    chargelist[iB_3_1] = charge[seqstring[iseq+1]]
    chargelist[iB_3_2] = charge[WC[seqstring[iseq+1]]]
    chargelist[iS_3_2] = charge["S"]
    return [coordlist[iS_5_1], coordlist[iB_5_1], coordlist[iB_5_2], coordlist[iS_5_2], coordlist[iP_1], coordlist[iP_2], coordlist[iS_3_1], coordlist[iB_3_1], coordlist[iB_3_2], coordlist[iS_3_2]]


def create_coordinates(seqstring, datafile):
	lseq = len(seqstring)
	seqstring = "X"+seqstring           ## Add leading "X" so that counting is easier
	natoms = lseq*6-2
	coordlist = np.zeros((natoms+1,3))
	typelist = np.zeros(natoms+1, dtype=int)
	chargelist = np.zeros(natoms+1)

	nbonds=2
	nangles=2
	ndihedrals=2

	file_writer=open(datafile,'w')

	#### Write Coordinates
	coord_mat_old = np.zeros((10,3))
	coord_mat_new = np.zeros((10,3))
	coord_ref = np.zeros((4, 3))
	coord_trj = np.zeros((4, 3))

	if lseq == 1:
	    iseq = 1
	    iS_5_1 = iseq*3-2
	    loc = copy.deepcopy(pos_dic[seqstring[1]+"A"])
	    z0 = loc["S_5_1"][2]
	    for loc_el in loc:
	        loc[loc_el][2] = loc[loc_el][2] - z0
	    iB_5_1 = iS_5_1 + 1
	    iB_5_2 = natoms + 1 - iS_5_1
	    iS_5_2 = natoms - iS_5_1
	    coordlist[iS_5_1] = loc["S_5_1"]
	    coordlist[iB_5_1] = loc["B_5_1"]
	    coordlist[iB_5_2] = loc["B_5_2"]
	    coordlist[iS_5_2] = loc["S_5_2"]
	    typelist[iS_5_1] = tipo["S"]
	    typelist[iB_5_1] = tipo[seqstring[iseq]]
	    typelist[iB_5_2] = tipo[WC[seqstring[iseq]]]
	    typelist[iS_5_2] = tipo["S"]
	else:
	    iseq = 1
	    loc = copy.deepcopy(pos_dic[seqstring[iseq:(iseq+2)]])
	    r0 = loc["S_5_1"]
	    #### Shift all coordinates so that first sugar has r=0
	    icoord = 0
	    for loc_el in loc:
	        loc[loc_el] = loc[loc_el] - r0
	        #coord_mat_old[icoord] = loc[loc_el]
	        icoord += 1
	    coord_mat_old = update_coordinates(loc, coordlist, typelist, chargelist, iseq, seqstring, 'y')
	    for iseq in range(2, lseq):
	        loc = copy.deepcopy(pos_dic[seqstring[iseq:(iseq+2)]])
	
	        ### create matrix with coordinates
	        icoord = 0
	        for loc_el in loc:
	            coord_mat_new[icoord] = loc[loc_el]
	            icoord += 1
	        
	        ### create ref and trj and align them
	        coord_ref = coord_mat_old[6:10]
	        coord_trj = coord_mat_new[0:4]
	        rmsd, U, trj_aligned = align(coord_ref, coord_trj)
	
	        ### use computed rotation matrix to align whole new step
	        centro_ref = np.mean(coord_ref, axis=0)
	        centro_trj = np.mean(coord_trj, axis=0)
	        coord_mat_new = coord_mat_new - centro_trj
	        coord_mat_new = np.dot(coord_mat_new, U)
	        coord_mat_new = coord_mat_new + centro_ref
	
	        ### create new dictionary
	        loc = {'S_5_1': coord_mat_new[0], 'B_5_1': coord_mat_new[1], 'B_5_2': coord_mat_new[2], 'S_5_2': coord_mat_new[3], 'P_1': coord_mat_new[4], 'P_2': coord_mat_new[5], 'S_3_1': coord_mat_new[6], 'B_3_1': coord_mat_new[7], 'B_3_2': coord_mat_new[8], 'S_3_2': coord_mat_new[9]}
	        
	        ### update coordinates
	        coord_mat_old = update_coordinates(loc, coordlist, typelist, chargelist, iseq, seqstring, 'n')
	
	
	    ### Make overall rotation such that line connecting first and last sugar is parallel to z axis
	    pos0 = coordlist[1,0:3]
	    pos1 = coordlist[int(0.5*natoms-1),0:3]
	    vec2 = pos1 - pos0
	    vec2 = vec2/np.linalg.norm(vec2)
	    vec1 = np.array([0,0,1])
	    rotvec = np.cross(vec2, vec1)
	    theta = acos(np.dot(vec2, vec1))
	    rotvec = rotvec/np.linalg.norm(rotvec)
	    for i in range(1,natoms+1):
	        vstart = coordlist[i,0:3] - coordlist[0,0:3]
	        pscal = np.dot(rotvec, vstart)
	        pvec = np.cross(rotvec, vstart)
	        vstop = vstart*np.cos(theta) + pvec*np.sin(theta) + rotvec*pscal*(1-np.cos(theta))
	        coordlist[i,0:3] = coordlist[0,0:3] + vstop
	
	for i in range(1,natoms+1):
	    print(str(i-1)+" "+"{:.4f}".format(coordlist[i,0])+" "+"{:.4f}".format(coordlist[i,1])+" "+"{:.4f}".format(coordlist[i,2]), file=file_writer)
	return coordlist, typelist, chargelist
