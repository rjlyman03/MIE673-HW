# Import Statements
import os
import matplotlib
import numpy as np
import math



def straight_beam_inertia(z,m):
    #z is an array of distances (m)
    #m is an array of linear mass densitys corresponding to z (kg/m)    
    M = np.trapezoid(m, z)
    S = np.trapezoid(m*z, z)
    J = np.trapezoid(m*(z**2), z)
    return M , S , J

def drive_train_inertia(n, J_h, J_g):
    J_dt = J_h + J_g * (n**2)
    return J_dt

def tower_ODE(m_t, k_t, zeta=0.33):
    omega = np.sqrt(k_t/m_t)
    return omega_n, c_t
