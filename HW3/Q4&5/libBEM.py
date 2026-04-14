import numpy as np
import pandas as pd
import math
#try:
from numpy import trapezoid
#except:
#    from numpy import trapz as trapezoid
import matplotlib
import matplotlib.pyplot as plt


def aero_coeffs(alpha, phi, fPolars):
    """
    Aerodynamic coefficients Cl, Cd, Cm, cn, and ct based on polar data and angle of attack.

    INPUTS:
     - alpha: array of angle of attack [deg or rad, consistent with fPolar ]
     - phi  : aray of flow angles [rad]
     - fPolars: list of polar interpolants or single polar interpolant such that 
                           Cl, Cd, Cm = fPolar(alpha) 
                           Cl, Cd, Cm = fPolars[0](alpha[0]) 
    OUTPUTS:
     - Cl, Cd, Cm, cn, ct: arrays of same shape as alpha 
    """
    if not hasattr(fPolars, '__len__'):
        fPolars = [fPolars] * len(alpha)
    Cl = np.zeros_like(alpha)
    Cd = np.zeros_like(alpha)
    Cm = np.zeros_like(alpha)

    #Cl, Cd, Cm = np.array([fPolar(alph) for fPolar, alph in zip(fPolars, alpha)]).T
    Cl, Cd, Cm = np.array([fPolar(alpha) for fPolar, alpha in zip(fPolars, alpha)]).T
    #Cl = 3.9*Cl
    cn1 = 3.25
    cn2 = 1.1
    cn1 = 1
    cn2 = 1
    #Cd = 1.1*Cd
    ct1 = 1
    ct2 = 1
    cn = Cl * cn1* np.cos(phi) + Cd * cn2 * np.sin(phi)
    ct = Cl * ct1* np.sin(phi) - Cd * ct2 * np.cos(phi)
    return Cl, Cd, Cm, cn, ct

def spanloads(pitch, Omega, U, r, uin, uit, chord, twist, fPolars, rho=1.225, B=3):
    """
    Computes spanwise aerodynamic loads.
    INPUTS:
     - pitch  (float): pitch angle of the blade [deg]
     - Omega  (float): rotational speed of the blade [rad/s]
     - U      (float): axial wind speed  [m/s]
     - r      (array): radial coordinates, measured from rotor center [m]
     - uin    (array): axial inductions, shape of r  [m/s]
     - uit    (array): axial inductions, shape of r  [m/s]
     - chord  (array): blade chord, shape of r  [m]
     - twist  (array): blade twist, shape of r  [deg]
     - fPolars(misc ): list of polar interpolants or single polar interpolant 
     - rho    (float): air density [kg/m^3]
     - B      (int  ): number of blades

    OUTPUTS:
     - df: dataframe with spanwise data, in particular, spanwise loads, fn, ft, mz
    """
    # --- Derived params
    rhub, R  = r[0], r[-1]

    # --- Flow variables
    #Using Textbook equations to solve for values:
    #we are assuming no induction, so a = a' = 0 => and uit = uin = 0
    Un   = U + uin
    Ut   = Omega*r + uit
    Vrel = np.sqrt((Un**2)+(Ut**2))
    phi  = np.arctan2(Un, Ut) #rad

    # --- Aerodynamic coefficients
    alpha_deg = phi*180/np.pi - (pitch + twist)
    alpha = np.deg2rad(alpha_deg)
    Cl, Cd, Cm, cn, ct = aero_coeffs(alpha_deg, phi, fPolars)
    #Cl, Cd, Cm, cn, ct = aero_coeffs(alpha, phi, fPolars)
    
    # --- Aerodynamic section loads: normal, tangential, and torsional moment
    # TODO 
    
    #since the units are given per meter, I am going to assume these should
    #be force/moment per unit span, distributed loads
    fD = 0.5*rho*(Vrel**2)*chord*Cd  #N/m
    fL = 0.5*rho*(Vrel**2)*chord*Cl  #N/m
    mz = 0.5*rho*(Vrel**2)*(chord**2)*Cm  #Nm/m 
   
    fn = fL*np.cos(phi) + fD*np.sin(phi) #N/m
    ft = fL*np.sin(phi) - fD*np.cos(phi) #N/m
    #fn = fL*np.cos(phi) + fD*np.sin(phi) #N/m
    #ft = fL*np.sin(phi) - fD*np.cos(phi) #N/m


    # --- Output
    data = [
        r, chord, Vrel, phi*180/np.pi, alpha_deg,
        Cl, Cd, cn, ct,  
        fn, ft, mz, 
        Un, Ut, uin, uit
    ]
    cols = [
        'r_[m]', 'chord_[m]', 'Vrel_[m/s]', 'phi_[deg]', 'alpha_[deg]',
        'Cl_[-]', 'Cd_[-]', 'cn_[-]', 'ct_[-]',
        'fn_[N/m]', 'ft_[N/m]', 'mz_[Nm/m]', 
        'Un', 'Ut', 'uin', 'uit'
    ]
    return pd.DataFrame(np.column_stack(data), columns=cols)

def intloads(r, fn, ft, mz, R, U, Omega, rho=1.225, B=3):
    """
    Computes integrated loads: thrust, torque, flapping moment and dimensionless coefficients CT, CQ, CP
    based on spanwise loads
    """
    r  = np.asarray(r)
    fn = np.asarray(fn)
    ft = np.asarray(ft)
    mz = np.asarray(mz)

    dT = B * fn 
    dQ = B * ft * r 
    dM = B * mz 
    
    T = trapezoid(dT, r)
    Q = trapezoid(dQ, r)
    M = trapezoid(dM, r)

    Mflap = trapezoid(r*fn, r) 
    Medge = trapezoid(r*ft, r)
    
    CT = T / (0.5 * rho * U**2 * np.pi * R**2)
    CQ = Q / (0.5 * rho * U**2 * np.pi * R**3)
    CP = CQ * Omega * R / U
    
    d={'T':T,'Q':Q, 'P':Q*Omega ,'M':M, 'CT':CT,'CQ':CQ, 'CP':CP, 'Mflap':Mflap, 'Medge':Medge}
    return d

def BEMqs(pitch, Omega, U, r, chord, twist, fPolars, rho=1.225, B=3, cone=0, a=None, ap=None):
    """
    Quasi-Steady Blade Element Momentum (BEM) calculation.
    
    INPUTS:
     - pitch  (float): pitch angle of the blade [deg]
     - Omega  (float): rotational speed of the blade [rad/s]
     - U      (float): axial wind speed  [m/s]
     - r      (array): radial coordinates, measured from rotor center [m]
     - chord  (array): blade chord, shape of r  [m]
     - twist  (array): blade twist, shape of r  [deg]
     - fPolars(misc ): list of polar interpolants or single polar interpolant 
     - rho    (float): air density [kg/m^3]
     - B      (int  ): number of blades
     - a      (array): initial guess for induction factors [-]
     - ap     (array): initial guess for induction factors [-]
    
    OUTPUTS:
     -  uin (array): Axial inflow velocity at each radial position
     -  uit (array): Tangential inflow velocity at each radial position
     -  dfB (dataframe): spanwise outputs
    """
    debug=True
    debug_r_section = 2
    a_list = []
    ap_list =[]
    iterations_a = []
    cCone  = np.cos(cone*np.pi/180.) #  = dr/dz (if no sweep)
    rPolar = r * cCone
    R = r[-1]
    # Initial guess for induction factors
    if a is None:
        a  = np.ones(len(r))*0.1
    if ap is None:
        ap = np.ones(len(r))*0.01

    nItMax=500
    aTol = 10**-6
    relaxation=0.5
    lambda_r = Omega * r / U

    for iterations in np.arange(nItMax):
        # --- Flow variables
        # Velocities and angles
        if debug:
            a_list.append(a[math.ceil((len(a)/debug_r_section))])
            #print(len(a)/debug_r_section)
            #print("How far along are we")
            ap_list.append(ap[math.ceil((len(a)/debug_r_section))])
            iterations_a.append(iterations)
            '''
        if min(ap)<0:
            print("a'<0")
        if max(a)<1:
            print("a<1")
            '''
        uin = -1*U*a #-2.5
        uit = Omega*r*ap
        Un   = U + uin
        Ut   = Omega*r + uit
        Vrel = np.sqrt((Un**2)+(Ut**2))
        phi  = np.arctan2(Un, Ut)# [rad]
        
        # --- Aerodynamic coefficients
        alpha_deg = phi*180/np.pi - (pitch + twist)
        #print(alpha_deg)
        alpha = np.deg2rad(alpha_deg)
        Cl, Cd, Cm, cn, ct = aero_coeffs(alpha_deg, phi, fPolars)
        #cn = cn+1.1  #delete this, hypothetical correction factor

        # --- Tip losse
        for i in range(len(phi)):
            if np.sin(phi[i])<0.01:
                phi[i] = 0.01
            
        F = (2/np.pi)*np.arccos(np.exp(-B*(R-r)/(2*r*np.sin(phi))))
        #print(F)
        F[F<=0]=0.5 # To avoid singularities

        # --- Induction factors with high thrust corrections
        a_last  = a
        ap_last = ap
    
        eps = 10**-8
        sigma = chord*B/(2*np.pi*r)  # NOTE: based on polar radial coordinate, rPolar
        a  = 1/(1 + (4*F*(np.sin(phi))**2)/(sigma*(cn + eps)))
        # for i in range(len(a)):
        #     if a[i] > 0.3: 
        #         fG = (6 - 3*a[i])/4
        #         Ct = cn[i]*sigma[i]*(Vrel[i]**2)/(U**2)
        #         a[i] = Ct/(4*F[i]*(1 - fG*a[i]))

        a = a*relaxation + (1 - relaxation)*a_last
        ap = (sigma*ct*Vrel**2)/(4*(1-a)*(U**2)*lambda_r)

        if debug and Omega==(11*U/R):
            print("Iteration: ", iterations, "  cn: ", cn[math.ceil((len(cn)/debug_r_section))])

        # --- Convergence check
        if (iterations > 3 and (np.mean(np.abs(a-a_last)) + np.mean(np.abs(ap-ap_last))) < aTol): 
            break
    if iterations == nItMax-1:
        print('Maximum iterations reached: Omega=%.2f V0=%.2f TSR=%.2f Pitch=%.1f' % (Omega, U, Omega*R/U, pitch))
    # --- Outputs
    
    if debug and Omega==(11*U/R):
    
        #a and a' convergence visual
        plt.figure()
        plt.plot(iterations_a,a_list,color='red',label='a')
        plt.plot(iterations_a,ap_list,color = 'blue', label="a'")
        plt.xlabel("iteration")
        plt.ylabel("a or a'")
        plt.title("induction factors for r = 90m, convergence")
        plt.legend()
        plt.show()    
        
    
    data = [
        r, alpha_deg, phi*180/np.pi, chord,
        Cl, Cd, cn, ct,  
        Vrel, Un, Ut, uin, uit,
        a, ap,
        F,
    ]
    cols = [
        'r_[m]', 'alpha_[deg]', 'phi_[deg]', 'chord_[m]',
        'Cl_[-]', 'Cd_[-]', 'cn_[-]', 'ct_[-]',
        'Vrel_[m/s]', 'Un_[m/s]', 'Ut_[m/s]', 'uin_[m/s]', 'uit_[m/s]',
        'a_[-]', 'ap_[-]',
        'F_[-]',
    ]
    dfB = pd.DataFrame(np.column_stack(data), columns=cols)

    return uin, uit, dfB



