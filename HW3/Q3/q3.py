## Import Statements
from scipy.interpolate import interp1d 
import sys
import numpy as np
import pandas as pd

B = 3 # Number of blades
lamda_design = 11
R = 185 # [m] radius
alpha_design = 10 # ax

## Section A
df = pd.read_csv('../data/ffa-w3-211_Re10M.csv') # load into dataframe
df.columns = ['alpha', 'C_l', 'C_d', 'C_m']

x = df['alpha'].values # shape inputs 
y = df[['C_l', 'C_d', 'C_m']].values # shape inputs
fPolar = interp1d(x, y, axis=0) # interpolant

C_l_design, C_d_design, C_m_design = fPolar(alpha_design)
print(f"{C_l_design=:.3f}, {C_d_design=:.3f}, and {C_m_design=:.3f}")

## Sections B & C
sys.path.extend(["./PlotPlanform.py"])
import PlotPlanform

## Section D
## PSUEDO CODE
##
''''
    Here we assume aero data is organized in an 2 dimensional array,
    one dimension being the distance r from the root and the other
    being the angles of attack of the aerofoil at that r. The other
    assumption is that chord length being as close to the optimal
    value is better than some other locally optimal solution given
    the constraints of a maximum chord length.
''' 

def optiaml_aero(Cd_data, Cl_data, alpha_data, r):
    idx = 0
    L_over_D = r*0
    Cl = r*0
    Cd = r*0 
    alpha = r*0
    for idx in range(len(r)):
        Cl_r_data = Cl_data[idx,:]
        Cd_r_data = Cd_data[idx,:]
        alpha_r_data = alpha_data[idx,:]
        L_over_Ds= Cl_r_data/Cd_r_data
        alpha_idx = np.argmax(L_over_Ds)
        Cl[idx] = Cl_r_data[alpha_idx]
        Cd[idx] = Cd_r_data[alpha_idx]
        alpha[idx] = alpha_r_data[alpha_idx]
    return Cl, Cd, alpha

def planform_variable_aero(r, R, Cl, Cd, alpha, root_c, max_c):
    chord, phi, a, ap = PlotPlanform.planform_ClCd(r, R,  TSR, Cl, Cd, B=3)
    chord[0] = root_c
    chord[chord>=max_c] = max_c
    return chord, phi, a, ap

