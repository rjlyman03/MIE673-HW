## Import Statements
from scipy.interpolate import interp1d 
import numpy as np
import pandas as pd

if __name__ == '__main__':
    # Code here is not executed in main
    print("testing...")
## Code
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

## Section D
## PSUEDO CODE
def optiaml_aero(Cd_data, Cl_data, alpha_data, r):
    idx = 0
    L_over_D = r*0
    Cl = r*0
    Cd = r*0 
    alpha = r*0
    for ri in r:
        Cl_r_data = Cl_data[idx,:]
        Cd_r_data = Cd_data[idx,:]
        alpha_r_data = alpha_data[idx,:]
        L_over_Ds= Cl_r_data/Cd_r_data
        alpha_idx = np.argmax(L_over_Ds)
        Cl[idx] = Cl_r_data[alpha_idx]
        Cd[idx] = Cd_r_data[alpha_idx]
        alpha[idx] = alpha_r_data[alpha_idx]
        idx += 1
    return Cl, Cd, alpha

def planform_variable_aero(Cl, Cd, alpha, root_c, max_c):
    
    return chord, phi, a, ap

