import numpy as np
from scipy.integrate import cumulative_trapezoid
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

import scipy.optimize as sciopt



def cumtrapz_flipped(v, z, initial=0):
    """ Flipped cumulative trapezoids (integrate from the end) """
    return np.flip( cumulative_trapezoid(np.flip(v), np.flip(z), initial=initial) )

def deflection(z, EI, p=None, Ftip=0, Mtip=0):
    """
    Compute deflection from a cantilever beam with stiffness EI(z) for a loading p(z)
    See lecture notes G.3

    INPUTS:
     - z: array of size(n) for position along the beam [m]
     - EI: array of size(n) for stiffness along the beam [Nm^2]
     - p: array of size(n) of normal load per unit span [N/m]
     - Ftip, Mtip: force and moment at tip
     
    OUTPUTS:
    - u, theta, kappa: deflection , slope and curvature
    - S, M: shear and bending moment
    """
    if p is None:
        p=z*0
    z = np.asarray(z)
    p = np.asarray(p)
    EI = np.asarray(EI)
    # --- Shear and moment
    S = cumtrapz_flipped(p, z, initial=0) - Ftip
    M = cumtrapz_flipped(S, z, initial=0) + Mtip
    # --- Curvature
    kappa = M / EI
    # --- Slope and deflection
    theta = cumulative_trapezoid(kappa, z, initial=0)
    u     = cumulative_trapezoid(theta, z, initial=0)
    return u, theta, kappa, S, M


def compute_modes(n_modes, z, EI, m, maxiter=500, abs_tol=0.01):
    """
    Use inverse iteration method to find mode shapes of a beam.
    See lecture notes G.4 and G.5
    INPUTS:
     - z: array of size(n) for position along the beam [m]
     - EI: array of size(n) for stiffness along the beam [Nm^2]
     - m: array of size(n) for mass per unit length along the beam [kg/m]
    OUTPUTS:
     - freqs: array of frequencies of size(n_modes) [Hz]
     - modes: array of mode shapes of size(n_modes, n) [-]
    """
    n = len(z)
    # --- Sanitization
    z = np.asarray(z)
    m = np.asarray(m)
    EI = np.asarray(EI)

    # --- Looping to find modes
    modes = np.zeros((n_modes, n))
    freqs = np.zeros(n_modes)
    for iMode in range(n_modes):
        pz = np.ones(n)
        omega2_prev = 0
        # --- Iteration loop
        for step in range(maxiter):
            uz , thetay, ky, Sy, Mz = deflection(z, EI, pz)
            # Orthogonalization against previous modes
            for iPrevMode in range(iMode):
                uz_p = modes[iPrevMode,:]
                const = trapezoid(uz_p * m * uz, z) / trapezoid(uz_p * m * uz_p, z)
                uz -=  const * uz_p
            # Normalization
            if uz[-1]==0:
                omega2, norm_factor = pz[-1] , 1   # Would unlikely happen
            else:
                omega2 = pz[-1] / (uz[-1] * m[-1])
                norm_factor = np.sqrt(uz[-1]**2)
            pz = omega2 * m * uz / norm_factor
            if abs(omega2 - omega2_prev) < abs_tol:
                #print(omega2, pz / (uz* m))
                break
            # Prepare for next iteration
            omega2_prev = omega2
        if step == maxiter - 1:
            raise Exception(f'Mode {iMode+1} did not converge')
        
        modes[iMode,:] = uz
        freqs[iMode]   = np.sqrt(omega2)/(2*np.pi)
        
    return freqs, modes

def generalized_MK(z, m, EI, phi, p_z=None, M_top=0, g=9.81):
    """
    Computes the generalized mass M (GM) and stiffnesses (GK_EI and GK_N) for a beam under
    axial load pz(z) and top mass Mtop for a given shape function Φ(z).

    INPUTS:
     - z: array of size(n) for position along the beam [m]
     - EI: array of size(n) for stiffness along the beam [Nm^2]
     - m: array of size(n) for mass per unit length along the beam [kg/m]
     - phi: array of size(n) for shape function along the beam [kg/m]
     - p_z: optional, array of size(n), axial distributed load along the beam [N/m]
     - M_top: optional, top mass above the beam [kg]
     - g: acceleration of gravity [m/s^2]
    OUTPUTS:
     - GM: generalized mass
     - GK_EI: baseline generalized stiffness 
     - GK_Np: generalized stiffness from axial loads p_z
     - GK_NM: generalized stiffness from top mass
    """
    # --- Sanitization
    z = np.asarray(z)
    m = np.asarray(m)
    EI = np.asarray(EI)
    phi = np.asarray(phi)
    # --- 1. Generalized mass
    GM = np.trapezoid(m * (phi**2), z)
    # --- 2. Generalized stiffnesses
    dphi = np.gradient(phi, z) # first derivative of phi wrt z
    d2phi = np.gradient(dphi, z) # second derivative of phi wrt z    

    GK_EI = np.trapezoid(EI * (d2phi**2), z)

    # G_KN
    GK_Np = 0
    if p_z is not None:
        p_z= np.asarray(p_z)
        N = -cumtrapz_flipped(p_z, z) # axial load
        GK_Np = trapezoid(N * dphi**2, z)
    GK_NM = -M_top * g * trapezoid(dphi**2, z)

    return GM, GK_EI, GK_Np, GK_NM




