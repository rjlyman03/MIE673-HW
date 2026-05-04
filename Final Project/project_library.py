import os
import matplotlib
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd
from sqlalchemy.engine import characteristics
import libBEM as bem
from scipy.interpolate import interp1d

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

def BEM(fPolars, pitch, Omega_rated, U_avg, rBEM, chord, twist, a, aprime, R_0):
    u_induced_normal, u_induced_tangential, DF1 = bem.BEMqs(pitch, Omega_rated, U_avg, rBEM, chord, twist, fPolars, rho=1.225, B=2, cone=0, a=a, ap=aprime)
    DF2 = bem.spanloads(pitch, Omega_rated, U_avg, rBEM, u_induced_normal, u_induced_tangential, chord, twist, fPolars, rho=1.225, B=2)
    f_normal, f_tangential, moment_aero_center = DF2["fn_[N/m]"].to_numpy(), DF2["ft_[N/m]"].to_numpy(), DF2["mz_[Nm/m]"].to_numpy()
    performanceDict = bem.intloads(rBEM, f_normal, f_tangential, moment_aero_center, R_0, U_avg, Omega_rated, rho=1.225, B=2)
    characteristicsDict = {
        "B": 2,
        "rho": 1.225,
        "U": U_avg,
        "omega": Omega_rated,
        "R": R_0,
        "lambda": (Omega_rated* R_0 / U_avg)
    } 
    scalarDict = characteristicsDict | performanceDict
    vectorDict = {
        "r": rBEM, 
        "c": chord, 
        "alpha": DF1["alpha_[deg]"].to_numpy(), 
        "beta": twist,
        "phi": DF1["phi_[deg]"].to_numpy(),
        "F": DF1["F_[-]"].to_numpy(), 
        "a": DF1["a_[-]"].to_numpy(), 
        "a_prime": DF1["ap_[-]"].to_numpy(), 
        "u_in": u_induced_normal, 
        "U_n": DF1["Un_[m/s]"].to_numpy(), 
        "u_it": u_induced_tangential, 
        "U_t": DF1["Ut_[m/s]"].to_numpy(), 
        "f_n": f_normal, 
        "f_t": f_tangential, 
        "m_aero":moment_aero_center
    }
    outDict = { "system": scalarDict, "spanwise": vectorDict }
    return outDict
