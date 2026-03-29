import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 13


def planform_nowakerot(r, R, TSR, Cl, B=3):
    """ 
    Chord and twist distribution for ideal rotor without wake rotation

    INPUTS:
     - r   (array): radial positions [m]tj
     - R   (float): rotor radius [m] 
     - TSR (float): design tip-speed ratio=Omega R/U0 [-]
     - Cl  (float): design lift coefficient [-]
     - B   (int  ): number of blades [-]

    OUTPUTS:
     - chord (array): chord [m] (shape of r)
     - phi   (array): twist [def] (shape of r)
     - a     (array): axial induction [-] (shape of r)
     - ap    (array): tangential induction [-] (shape of r)
    """
    r = np.asarray(r)
    lambda_r = TSR*r/R   
    phi   = np.arctan(2/(3*lambda_r))
    a     = np.ones_like(r)/3
    chord = 8*np.pi*r*np.sin(phi)/(3*B*Cl*lambda_r)
    ap    = r*0

    phi *= 180/np.pi            # [deg]
    return chord, phi, a, ap


def planform_wakerot(r, R, TSR, Cl, B=3):
    """ 
    Chord and twist distribution for ideal rotor with wake rotation

    INPUTS:
     - r   (array): radial positions [m]
     - R   (float): rotor radius [m] 
     - TSR (float): design tip-speed ratio=Omega R/U0 [-]
     - Cl  (float): design lift coefficient [-]
     - B   (int  ): number of blades [-]

    OUTPUTS:
     - chord (array): chord [m] (shape of r)
     - phi   (array): twist [def] (shape of r)
     - a     (array): axial induction [-] (shape of r)
     - ap    (array): tangential induction [-] (shape of r)
    """
    # TODO
    lambda_r = TSR*r/R
    phi = (2/3)*np.arctan(1/lambda_r)
    print(phi)
    chord = 8*np.pi*r*(1 - np.cos(phi))/(B*Cl) 
    sigma = B*chord/(2*np.pi*r)
    a     = 1/(1 + (4*np.sin(phi)**2)/(sigma*Cl*np.cos(phi)))
    ap    = (1 - 3*a)/(4*a - 1)
    phi *= 180/np.pi            # [deg]
    return chord, phi, a, ap


def planform_ClCd(r, R,  TSR, Cl, Cd, B=3):
    """ 
    For Grad/Groups
    """
    lambda_r = TSR*r/R

	# Axial induction
    bValid = lambda_r>0
    a = np.ones_like(r) * 1/4 # Initialization, with value for lambda_r ==0
    a[bValid] = 1/2 * (1 - np.sqrt(1 + lambda_r[bValid]**2) * np.sin((1/3)*np.arctan(1/lambda_r[bValid])))
 	
 	# Tangential induction
    ap = np.full_like(lambda_r, np.nan)  # Initialization, with value for lambda_r ==0
    ap[bValid] = (1 - 3*a[bValid])/(4*a[bValid] - 1)

 	# Flow angle and chord accounting for tip loss and drag
    phi   = np.arctan((1 - a)/((1 + ap)*lambda_r))
    Cn    = np.cos(phi)*Cl + np.sin(phi)*Cd # TODO
    F     = (2/np.pi)*np.arccos(np.exp(-B*(R-r)/(2*r*np.sin(phi))))
    chord = a*(8*np.pi*r*F*np.sin(phi)**2)/((1-a)*B*Cn) # TODO
    phi *= 180/np.pi            # [deg]
    return chord, phi, a, ap 


''' if __name__ == '__main__': '''
testing=1
if (testing==1):
    B            = 3   # Number of blades
    TSR          = 11  # Design tip speed ratio [-]
    alpha_design = 10  # Design angle of attack [deg]
    R            = 185 # Total tip radius [m]
    r_hub        = 5   # Hub radius [m]
    r   = np.linspace(r_hub, R, 100)

    Cl    = 1.51 # based on polar
    Cd    = 0.14 # based on polar

    # No Wake Rot
    chord_nw, phi_nw, a_nw, ap_nw = planform_nowakerot(r, R, TSR, Cl, B)
    twist_nw = phi_nw-alpha_design
    # Wake Rot
    chord_wr, phi_wr, a_wr, ap_wr = planform_wakerot(r, R, TSR, Cl, B)
    twist_wr = phi_wr-alpha_design

    # Wake Rot and Tip Loss
    chord_ld, phi_ld, a_ld, ap_ld = planform_ClCd(r, R, TSR, Cl, Cd*10, B)
    twist_ld = phi_ld-alpha_design

    # --- Plot
    fig,axes = plt.subplots(2, 2, sharex=True, figsize=(12.8,4.8))
    fig.subplots_adjust(left=0.09, right=0.99, top=0.97, bottom=0.12, hspace=0.30, wspace=0.25)
    
    # Chord
    ax = axes[0,0]
    ax.plot(r/R, chord_nw/np.max(r) ,'-'  , label='Without wake rotation')
    ax.plot(r/R, chord_wr/np.max(r) ,'--' , label='With wake rotation')
    ax.plot(r/R, chord_ld/np.max(r) ,'.'  , label='With wake rotation + Drag + Tiploss')
    ax.set_ylabel(r'Blade chord, $c/R$ [-]       ')
    ax.legend()
    
    # Twist
    ax = axes[0,1]
    ax.plot(r/R, twist_nw, '-'  , label='Without wake rotation')
    ax.plot(r/R, twist_wr, '--' , label='With wake rotation')
    ax.plot(r/R, twist_ld ,'.'  , label='With wake rotation + Drag + Tiploss')
    ax.set_ylabel('Blade twist angle [deg]')
    
    # Axial induction
    ax = axes[1,0]
    ax.plot(r/R, a_nw, '-'  , label='Without wake rotation')
    ax.plot(r/R, a_wr, '--' , label='With wake rotation')
    ax.plot(r/R, a_ld ,'.'  , label='With wake rotation + Drag + Tiploss')
    ax.set_xlabel(r'Blade radius, $r/R$ [-]')
    ax.set_ylabel('Axial induction [-]')
    
    # Tangential induction
    ax = axes[1,1]
    ax.plot(r/R, ap_nw, '-'  , label='Without wake rotation')
    ax.plot(r/R, ap_wr, '--' , label='With wake rotation')
    ax.plot(r/R, ap_ld ,'.'  , label='With wake rotation + Drag + Tiploss')
    ax.set_xlabel(r'Blade radius, $r/R$ [-]')
    ax.set_ylabel('Tangential induction [-]')

    # --- add ref turbines
    #from openfast_toolbox.io import FASTInputFile
    #import os
    #scriptDir = os.path.dirname(__file__)
    #stys = ['-','--',':','.-','o']
    #vP_rated = [15, 22]
    #for i, pr in enumerate(vP_rated): 
    #    folder='{:04.1f}MW'.format(pr)
    #    df = FASTInputFile(os.path.join(scriptDir, '../../data-ref/'+folder+'/AeroDyn_blade.dat')).toDataFrame()
    #    # --- Plot
    #    r1 = df['BlSpn_[m]']
    #    c1 = df['BlChord_[m]']
    #    t1 = df['BlTwist_[deg]']
    #    axes[0,0].plot(r1/np.max(r1), c1/np.max(r1), stys[i], label='{:.2f}MW'.format(pr))
    #    axes[0,1].plot(r1/np.max(r1), t1, stys[i], label='{:.2f}MW'.format(pr))
    #axes[0,0].legend()
    #axes[0,1].legend()
    #fig.savefig('planforms.png')

    plt.show()
