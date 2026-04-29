import os
import matplotlib
import matplotlib.pyplot as plt
# from openfast_toolbox.io import FASTInputFile
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
