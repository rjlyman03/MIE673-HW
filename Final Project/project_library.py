import os
import matplotlib
import matplotlib.pyplot as plt
# from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd
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

def BEM(DF0,pitch,Omega_rated,U_avg,rBEM,chord,twist,a,aprime,R_0):
    A = np.array(DF0)
    alpha_polar = A[:,0]
    coefs_polar = A[:,1:]
    #print(alpha_polar.shape, coefs_polar.shape)
    fPolars = interp1d(alpha_polar, coefs_polar, axis=0)
    u_induced_normal, u_induced_tangential, DF1 = bem.BEMqs(pitch, Omega_rated, U_avg, rBEM, chord, twist, fPolars, rho=1.225, B=2, cone=0, a=a, ap=aprime)
    DF2 = bem.spanloads(pitch, Omega_rated, U_avg, rBEM, u_induced_normal, u_induced_tangential, chord, twist, fPolars, rho=1.225, B=2)
    f_normal, f_tangential, moment_aero_center = DF2["fn_[N/m]"].to_numpy(), DF2["ft_[N/m]"].to_numpy(), DF2["mz_[Nm/m]"].to_numpy()
    DICT = bem.intloads(rBEM, f_normal, f_tangential, moment_aero_center, R_0, U_avg, Omega_rated, rho=1.225, B=2)
    return DICT,DF1,DF2,f_normal,f_tangential
