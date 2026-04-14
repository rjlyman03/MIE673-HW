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
p.L_t =      149.39 # m
p.M_t =      1.578 * 10**6 # kg
p.S_t =      9.714 * 10**7 # kg.m
p.J_t =      8.645 * 10**9 # kg.m^2
p.mtt =      TODO # kg    Generalized mass
p.k_t =      TODO # N/m   Generalized stiffness (for verif: ~3e+06)
p.c_t =      TODO # N/m.s Generalized damping
# --- Blade 
p.L_b =   137.80 # m
p.M_b = 82427.56 # kg
p.S_b = 3.07e+06 # kg.m
p.J_b = 2.08e+08 # kg.m^2
p.mbb =  2910.13 # kg    Generalized mass
p.k_b = 18180.05 # N/m   Generalized stiffness
p.c_b =    72.74 # N/m.s Generalized damping
# --- Hub/Shaft 
p.M_s = 1.20e+05 # kg
p.J_s = 1.88e+06 # kg.m^2
# --- Drivetrain (Hub+Gen LSS)
p.J_DT = TODO 
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
# ED     = FASTInputFile('../../data-ref/22.0MW/ElastoDyn.dat')
# ed_bld = FASTInputFile('../../data-ref/22.0MW/ElastoDyn_blade.dat').toDataFrame()
# ad_bld = FASTInputFile('../../data-ref/22.0MW/AeroDyn_blade.dat').toDataFrame()
# polar  = FASTInputFile('../../data-ref/22.0MW/Airfoils/IEA-22-280-RWT_AeroDyn15_Polar_50.dat').toDataFrame().values[:,:4]
# print(ed_bld.keys()) 
# print(ad_bld.keys())
# # Blade shape function
# p.z_b   = np.asarray(ed_bld['BlFract_[-]']) * p.L_b + ED['HubRad']
# p.phi_b = np.asarray(ed_bld['ShapeFlap1_[-]'])
# # Blade Aerodynamics
# p.r_full = np.asarray(ad_bld['BlSpn_[m]']) + ED['HubRad'] # We are lucky it's the same as p.z_b
# p.twist  = np.asarray(ad_bld['BlTwist_[deg]'])
# p.chord  = np.asarray(ad_bld['BlChord_[m]'])
# p.fPolar = interp1d(polar[:,0], polar[:,1:], axis=0) 
# np.testing.assert_array_almost_equal(p.z_b, p.r_full)






# --------------------------------------------------------------------------------
# --- 4.5b State-space model
# --------------------------------------------------------------------------------
def statespace(t, x, p, fu=None, calcOutput=False):
    """
    First-order state-space model for the 3DOF wind turbine
    
    INPUTS:
     - x: 1d array [x_n, psi, x_b, xdot_n, psi_dot, xdot_b]
     - fu: function handle u(t)
     - p: name space with parameters
     - calcOutput: if True, return additional outputs
    
    OUTPUT:
     - xdot, 1d array:  [xdot_n, psi_dot, xdot_b, xddot_n, psi_ddot, xddot_b]
     - out: optional, pandas series
    """
    # --- Outputs
    xdot = np.zeros_like(x)
    out = {'t':t, 'x_n':x[0], 'psi':x[1], 'x_b': x[2], 'xd_n':x[3], 'psid':x[4], 'xd_b':x[5]}

    # --- Aliases if needed
    x_n = x[0]
    psi = x[1]

    # --- Generator torque
    Qg = 0
    out['Qg'] = Qg

    # --- Aerodynamic calculations
    if fu is not None:
        U0 = fu(t)
        out['U0'] = U0 # Store additional time outputs

    # --- Simple transfer of derivatives
    xdot[0] = x[3] # xdot_n = xdot_n
    # TODO TODO TODO

    # --- Main ODE:
    xdot[3] = -1/p.mnt * p.k_t * x_n # xddot_n (Dummy model)
    # TODO TODO TODO

    if calcOutput:
        out.update({'xdd_n':xdot[3], 'psidd':xdot[4], 'xdd_b':x[5]})
        out = pd.Series(out)
        return out
    else:
        return xdot

def calcOutput(res, p, fu):
    """ Calculate output based on knowledge of the state """
    rows = []
    for ti, xi in zip(res.t, res.y.T):
        out_i = statespace(ti, xi, p, fu=fu, calcOutput=True)
        rows.append(out_i)
    dfout = pd.DataFrame(rows)
    return dfout


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





# ---
plt.show()
