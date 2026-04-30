import os
import matplotlib
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd

#wind
import numpy as np

#FIND R
def findR(P, rho, U, C_p):
    R = np.sqrt(P/(0.5*np.pi*rho*C_p*U**3))
    return R

def readPolarDatFile(filename):
    df = pd.read_table(filename, sep="\s+", names=["Alpha", "Cl", "Cd", "Cm"], skiprows=54, nrows=200)
    return df

def modifyBlade(r, chord, twist):
    fileToModify = "./25.0MW/AeroDyn_blade.dat"
    bld = FASTInputFile(fileToModify)

    nr = len(bld['BldAeroNodes'])
    print('Number of radial positions (SHOULD MATCH n) :', nr)

    bld['BldAeroNodes'] = np.zeros((nr, 7))

    #  0:   BlSpn    1: BlCrvAC ->0   2:BlSwpAC->0 3: BlCrvAng->0    4->BlTwist  5:BlChord          BlAFID
    bld['BldAeroNodes'][:,0] =  r   # r
    bld['BldAeroNodes'][:,1] =  r*0 #
    bld['BldAeroNodes'][:,2] =  r*0 #
    bld['BldAeroNodes'][:,3] =  r*0 #
    bld['BldAeroNodes'][:,4] = twist # twsit
    bld['BldAeroNodes'][:,5] = chord # chord
    bld['BldAeroNodes'][:,6] =  1 # ID

    bld.write(fileToModify)

    return # modify's .dat file

