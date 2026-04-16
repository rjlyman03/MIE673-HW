import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import types
from openfast_toolbox.io.fast_input_file import FASTInputFile
from scipy.interpolate import interp1d

TODO = 0
p = types.SimpleNamespace()
# ------------------------------------------------------------------
# --- 4.1 Main numerical parameters computed in homework
# ------------------------------------------------------------------
# --- Tower    
p.L_t =      149.39 # m                    #can add more sig figs
p.M_t =      1.578 * 10**6 # kg            #can add more sig figs
p.S_t =      9.714 * 10**7 # kg.m          #can add more sig figs
p.J_t =      8.645 * 10**9 # kg.m^2        #can add more sig figs
p.mtt =      2.390663 * 10**5  # kg    Generalized mass
p.k_t =      2.9265 * 10**6 # N/m   Generalized stiffness (for verif: ~3e+06)
p.c_t =      2.760243 * 10**5 # N/m.s Generalized damping
# --- Blade 
p.L_b =   137.80 # m
p.M_b = 82427.56 # kg
p.S_b = 3.07e+06 # kg.m
p.J_b = 2.08e+08 # kg.m^2
p.mbb =  2910.13 # kg    Generalized mass
p.mbt= 7575.00 # kg    
p.k_b = 18180.05 # N/m   Generalized stiffness
p.c_b =    72.74 # N/m.s Generalized damping
# --- Hub/Shaft 
p.M_s = 1.20e+05 # kg
p.J_s = 1.88e+06 # kg.m^2
# --- Drivetrain (Hub+Gen LSS)
p.J_DT = 5.000096e+06 #kg m^2 
# --- Nacelle
p.M_nac = 8.50e+05        # kg, Nacelle + Yaw
p.M_top = 1.05e+06        # kg, Nacelle + Yaw + Hub + Blade
p.mnt   = p.mtt + p.M_nac # kg, Nacelle + Generalized Tower 
# --- Environment
p.g = 9.81

# --- Derived parameters
zeta_b = 0.005 # [-]
f_b = 0.40     # [Hz]
T_b = 1/(f_b*np.sqrt(1-zeta_b**2))
print('Damped blade period {:.2f}'.format(T_b))


# --------------------------------------------------------------------------------}
# --- 4.5c Parameters needed for aerodynamics calculations
# --------------------------------------------------------------------------------{
ED     = FASTInputFile('../../mie673-data/22.0MW/ElastoDyn.dat')
ed_bld = FASTInputFile('../../mie673-data/22.0MW/ElastoDyn_blade.dat').toDataFrame()
ad_bld = FASTInputFile('../../mie673-data/22.0MW/AeroDyn_blade.dat').toDataFrame()
polar  = FASTInputFile('../../mie673-data/22.0MW/Airfoils/IEA-22-280-RWT_AeroDyn15_Polar_50.dat').toDataFrame().values[:,:4]
print(ed_bld.keys()) 
print(ad_bld.keys())
# Blade shape function
p.z_b   = np.asarray(ed_bld['BlFract_[-]']) * p.L_b + ED['HubRad']
p.phi_b = np.asarray(ed_bld['ShapeFlap1_[-]'])
# Blade Aerodynamics
p.r_full = np.asarray(ad_bld['BlSpn_[m]']) + ED['HubRad'] # We are lucky it's the same as p.z_b
p.twist  = np.asarray(ad_bld['BlTwist_[deg]'])
p.chord  = np.asarray(ad_bld['BlChord_[m]'])
p.fPolar = interp1d(polar[:,0], polar[:,1:], axis=0) 
np.testing.assert_array_almost_equal(p.z_b, p.r_full)





# --------------------------------------------------------------------------------
# --- 4.5b State-space model
# --------------------------------------------------------------------------------
x0 = [1, 0, 0, 0, 0, 0.1047] # initial conditions
def statespace(t, x, p):
    """
    First-order state-space model for the 3DOF wind turbine
    
    INPUTS:
     - x: 1d array [x_n, x_b, psi, dx_n, dx_b, dpsi]
     - p: name space with parameters
     - calcOutput: if True, return additional outputs
    
    OUTPUT:
     - xdot, 1d array:  [dx, dx_b, dpsi, ddx, ddx_b, ddpsi]
    """
    # --- Outputs
    x = np.asarray(x, 'C')
    xdot = np.zeros_like(x) 

    # --- Aliases for forces
    x_n = x[0]
    x_b = x[1]
    psi = x[2]
    dx_n = x[3]
    dx_b = x[4]
    dpsi = x[6]
    
    # -- Matrices
    M = np.array([[p.mnt, p.mbt, 0], [p.mbt, p.mbb, 0],[0, 0, p.J_b + p.J_DT ]])
    K = np.array([[p.k_t, 0, 0],[0, p.k_b, 0],[0, 0, 0]])
    A = np.block([np.zeros(3), np.eye(3), np.zeros(3), np.invert(M)*(-1*K)])
    
    # BEM Code + Generalized forces

    # --- Aerodynamic calculations
    U0 = 10
    lambda0 = 11
    P0 = 15.6 * 10**6
    Omega0 = 0.77
    Qg = (P0/Omego0) - (P0/Omego0)*(dpsi - Omega0)
    dummy1, dummy2, df = BEMqs(0, dpsi, U0, p.z_b, p.chord, p.twist, p.fPolar)
    cn = df["c_n"]
    ct = df["c_t"]
    R = 0.5*(p.z_b[:-2] + p.z_b[1:])
    fn = 0.5*1.225*(U0**2)*p.chord*cn
    FN = np.trapezoid(fn, R)
    FNB = np.trapezoid(fn*0.5*(p.phi_b[:-2] - p.phi_b[1:]), R)
    ft = 0.5*1.225*(U0**2)*p.chord*ct
    QA = np.trapezoid(ft*R, R)
    
    INT0 = np.trapezoid(p.phi_b, p.z_b)
    INT1 = p.S_b*np.trapezoid(np.gradient(p.phi_b, p.z_b), p.z_b) # wrong but close to write without p.m_prime
    INT2 = np.trapezoid(0.5*(p.phi_b[:-2] - p.phi_b[1:])*fn, R)
    INT3 = 0 # this one is really scuffed as it has m_prime*r*phi_b*phi_b_prime*
    Q1 = -p.c_t*dx_n - p.c_b*dx_b - FN + (dpsi**2) * dx_b * INT1
    Q2 = (-p.c_t*dx_n - p.c_b*dx_b)*INT0 + INT2 + INT3
    Q3 = -Qg + p.S_b*9.81*np.sin(psi) + QA
    B = np.invert(M)*np.array([0, 0, 0, Q1, Q2, Q3])

    # --- Simple transfer of derivatives
    xdot = A*x + B
return xdot

# --------------------------------------------------------------------------------
# --- 4.5b Numerical integration
# --------------------------------------------------------------------------------
x0 = np.array([0.5, 0, 0, 0, 0, 0]) # Initial conditions
fu = lambda t: 0*t
t = np.linspace(0, 60, 500)
res = solve_ivp(statespace, (t[0], t[-1]), x0, t_eval=t, args=(p, fu) )

# --- Plot
fig, axes = plt.subplots(5, 1, sharex=True, figsize=(12.8,6.0))
fig.subplots_adjust(left=0.07, right=0.98, top=0.98, bottom=0.08, hspace=0.16, wspace=0.20)
res.y[1,:] = np.mod(res.y[1,:], 2*np.pi)
dfOut = calcOutput(res, p, fu)

axes[0].plot(res.t, res.y[0]          , '-'); axes[0].set_ylabel(r'$x_n$ [m]')
axes[1].plot(res.t, res.y[1]*180/np.pi, '-'); axes[1].set_ylabel(r'$\psi$ [deg]')
axes[2].plot(res.t, res.y[2]          , '-'); axes[2].set_ylabel(r'$x_b$ [m]')
axes[3].plot(res.t, res.y[4]*30/np.pi , '-'); axes[3].set_ylabel(r'$\Omega$ [RPM]')
axes[4].plot(dfOut['t'], dfOut['Qg']  , '-'); axes[4].set_ylabel(r'$Q_g$ [Nm]')
axes[4].set_xlabel('Time [s]')


# --------------------------------------------------------------------------------}
# --- 4.5c 
# --------------------------------------------------------------------------------{

p.eff = 1
p.opt_psi_dot = (11*10)/p.r_full # omega
U0 = 10 # m/s, wind speed

p.opt_P = 15.8*10**-6 # MW, optimal power


p.Qg = ((1/p.eff)*(p.opt_P/p.opt_psi_dot) # - (1/p.eff)*(p.opt_P/(p.opt_psi_dot**2))) * (p.psi_dot - p.opt_psi_dot) - Taylor series explansion 

plt.show()

