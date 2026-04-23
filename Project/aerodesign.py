""" 
Functions for aerodynamic design

References:
   [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer
   [2] Manwel et al (2024) Wind Energy Explained, 3rd edition, Wiley
"""

import numpy as np

def planform_ClCd(r, R,  TSR_design, Cl_design, Cd_design, B=3):
    """ 
    Chord, and flow angle, assuming:
      - a, a' from Optimal AD theory
      - constant Cl, and Cd along the blade span

    NOTE:  a and aprime will be the same as optimal planform with wake rotation
           They are not changed by Cd or F because they are "actuator disc" inductions
    INPUTS:
      - r: radial position, array
      - R: rotor radius, scalar
      - TSR: tip speed ratio, scalar
      - Cl: Lift coefficient, scalar
      - Cd: Drag coefficient, scalar
      - B: Number of blades, scalar
    """

    lambda_r =r/R * TSR_design
    # --- Optimal AD results
    # Axial induction factor
    phi   = 2./3. * np.arctan(1./(lambda_r))
    #a = 1/(1+ np.sin(phi)**2/(  (1-np.cos(phi))*np.cos(phi)   ) )
    ## Tangential induction based on orthogonality of induced velocity for an AD
    #aprime  = (1 - 3 * a) / (4 * a - 1)                     # Only true for optimal rotor
    a, aprime = ADMTO_inductions(lambda_r)
    # --- General results
    cn     = Cl_design * np.cos(phi) + Cd_design * np.sin(phi)
    F      = (2/np.pi) * np.arccos(np.exp(-B / 2 * (R - r) / (r * np.sin(phi))))
    chord  = a/(1-a) * (8 * np.pi * F * r * (np.sin(phi))**2) / (B*cn)
    phi *= 180/np.pi            # [deg]

    return chord, phi,  a, aprime 

def ADMTO_inductions(lambda_, *args, **kwargs):
    """
    Actuator disc momentum theory *optimal* (ADMTO) induction factors
    See e.g. [1] Section 9.5.4
    """
    lambda_ = np.asarray(lambda_)
    a       = np.zeros_like(lambda_)
    ap      = np.zeros_like(lambda_)
    bValid = lambda_>0
    bZero  = lambda_==0
    a[bValid] = 1/2 * (1 - np.sqrt(1+lambda_[bValid]**2) * np.sin ( 1/3*np.arctan(1/lambda_[bValid]) )  ) # [1] Eq. 9.86
    #phi   = 2./3. * np.arctan(1./(lambda_[bValid]))
    #a[bValid] = 1/(1+ np.sin(phi)**2/(  (1-np.cos(phi))*np.cos(phi)   ) )
    # Limits:
    a[bZero] = 1/4
    ap[bValid] = (1-3*a[bValid])/(4*a[bValid]-1) # [1] Eq. 9.84, only true for optimal
    #ap= 1/2*(-1+np.sqrt(1+4/lambda_r**2 * a * (1-a)))    # Always true 
    ap[bZero]  = np.nan

    return a, ap

def CP_Lambda_AD(lambda_, method='analytical'):
    """
    Optimal power coefficient for a rotating actuator disc.
    See e.g. [2] Chap 3 
    INPUTS:
      - lambda_: tip speed ratio, array
    OUTPUT:
      - CP: power coefficient for each lambda_ values, array
    """
    lambda_ = np.asarray(lambda_)
    CP      = np.zeros_like(lambda_)
    a, ap = ADMTO_inductions(lambda_)

    def dCPmax(x): # [2] Chapter 3, eq 43
        return  (64/5*x**5 + 72*x**4 + 124*x**3 + 38*x**2 - 63*x-12*np.log(x) - 4/x)
    for i, (l,a2) in enumerate(zip(lambda_,a)):
        if l==0:
            CP[i]==0
        else:
            CP[i] = 8/(729*l**2) * ( dCPmax(x=1/4) - dCPmax(x=1-3*a2)  )

    return CP


def CP_lambda_B_ClCd_Wilson(lambda_, B, ClCd):
    """ 
    Optimial power coefficient for a rotor with B blades and a constant Cl over Cd
    Approximate formula from Wilson and Lissaman (1976)

    INPUTS: 
      - lambda_: tip speed ratio, array
      - B: number of blade, scalar
      - ClCd: lift over drag ration, scalar
    OUTPUT:
      - CP: power coefficient for each lambda_ values, array
    """
    lambda_ = np.asarray(lambda_)
    CP = np.zeros_like(lambda_)
    bValid = np.logical_and( lambda_>=4 , lambda_<=20)
    TSR = lambda_[bValid]
    T1 = TSR + (1.32+((TSR-8)/20)**2)/ B**(2/3)
    if ClCd>1e8:
        T2=0
    else:
        T2 = 0.57 * TSR**2 / (ClCd* (TSR+1/(2*B)) )
    CP[bValid] = 16/27 * TSR /T1 - T2
    CP[~bValid] = np.nan
    return CP
