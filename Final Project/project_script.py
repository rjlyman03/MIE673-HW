import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd
import aerodesign as ad
import libBEM as bem
import project_library as lib

#wind
import numpy as np

months = [
    "January","February","March","April","May","June",
    "July","August","September","October","November","December","Annual"
]

# Normal Monthly Average Speed [MPH]
avg_speed_mph = np.array([
    45.6, 44.6, 39.8, 35.6, 29.6, 26.8,
    25.5, 24.5, 27.6, 35.5, 39.4, 44.0, 34.9
], dtype=float)
mph2mps = 0.447
U_avg = mph2mps*np.mean(avg_speed_mph)
print(f"Yearly Average Wind Speed = {U_avg:,.4} meters per second " )

#Assumed Params
B = 2 # number of blades
alpha_design = 8.8 # degrees
C_l_design = 1.4
C_d_design = 0.012
Cl_by_Cd = C_l_design/C_d_design
C_P_rated_0 = 0.48
P = 25 * 10**6 # Watts
rho = 1.225 # wrong but it's fine
R_0 = lib.findR(P, rho, U_avg, C_P_rated_0)
print(f"Radius from Rated Power = {R_0:,.4} m")

# Do wind step and then turbSim!!!
test_lambdas = np.arange(4,20,0.1) # data has been fit for this range of TSR
C_P = ad.CP_lambda_B_ClCd_Wilson(test_lambdas, B, Cl_by_Cd)
index_lambda = np.argsort(np.abs(C_P - C_P_rated_0)) # indices of TSR, sorted from lowest to highest difference between TSR and the design TSR
lambda_design = test_lambdas[index_lambda[0]] # second index is the smaller TSR
print(f"Designed TSR = {lambda_design:,.4}")
Omega_rated = (lambda_design*U_avg)/R_0

n = 20 # number of elements
hubRad = 5 # a little bigger than the 22 MW
rBEM = np.linspace(hubRad, R_0, n + 1) # radial distances between elements
bladeFrac = rBEM/R_0 # length along the blade as a fraction of total blade length
chord, phi,  a, aprime = ad.planform_ClCd(rBEM, R_0, lambda_design, C_l_design, C_d_design, B)
twist = phi - alpha_design
pitch = 0*twist

#plot actuator disk w/ wake:
fig, axs = plt.subplots(4, sharex=True)
fig.suptitle("Optimal Actuator Disk Model with Wake Rotation")
fig.supxlabel("Blade Fraction")
ylabels = ["chord length (m)","flow angle (deg)","axial induction", "tangential induction"]
yvalues = [chord, phi, a, aprime]
idx=0
for ax in axs:
    ax.plot(bladeFrac, yvalues[idx])
    ax.set(ylabel=ylabels[idx])
    ax.grid(True)
    ax.label_outer()
    idx+=1
plt.legend()
plt.show()

# DO A BEM
DF0 = lib.readPolarDatFile("./22.0MW/Airfoils/IEA-22-280-RWT_AeroDyn15_Polar_56.dat")

"""
A = np.array(DF0)
alpha_polar = A[:,0]
coefs_polar = A[:,1:] 
print(alpha_polar.shape, coefs_polar.shape)
fPolars = interp1d(alpha_polar, coefs_polar, axis=0)
u_induced_normal, u_induced_tangential, DF1 = bem.BEMqs(pitch, Omega_rated, U_avg, rBEM, chord, twist, fPolars, rho=1.225, B=2, cone=0, a=a, ap=aprime)
DF2 = bem.spanloads(pitch, Omega_rated, U_avg, rBEM, u_induced_normal, u_induced_tangential, chord, twist, fPolars, rho=1.225, B=2)
f_normal, f_tangential, moment_aero_center = DF2["fn_[N/m]"].to_numpy(), DF2["ft_[N/m]"].to_numpy(), DF2["mz_[Nm/m]"].to_numpy()
DICT = bem.intloads(rBEM, f_normal, f_tangential, moment_aero_center, R_0, U_avg, Omega_rated, rho=1.225, B=2)
"""
DICT,DF1,DF2,f_normal,f_tangential = lib.BEM(DF0,pitch,Omega_rated,U_avg,rBEM,chord,twist,a,aprime,R_0)

print(DICT["CP"])
#DICT["CT"]
#DICT["CQ"]

#plot quasi-steady BEM
fig, axs = plt.subplots(3, sharex=True)
fig.suptitle("BEM loads and Angle of Attack")
fig.supxlabel("Blade Fraction")
ylabels = ["f_normal (N/m)","f_tangential (N/m)","AoA"]
yvalues = [f_normal, f_tangential, DF1["alpha_[deg]"]]
idx=0
for ax in axs:
    ax.plot(bladeFrac, yvalues[idx])
    ax.set(ylabel=ylabels[idx])
    ax.grid(True)
    ax.label_outer()
    idx+=1
plt.legend()
plt.show()

# MAKE THE ABOVE A FUNTION CALL IN THE LIBRARY

# Iterated close to Omega_rated
percent_away = 20
n_iterations = 5
Omega_upper = Omega_rated + (percent_away/100)*Omega_rated
Omega_lower = Omega_rated - (percent_away/100)*Omega_rated
for omega_i in np.linspace(Omega_lower, Omega_upper, n_iterations):
    print(omega_i)

