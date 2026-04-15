# Import Statements
import os
import matplotlib
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd
from libbeam import deflection,compute_modes,generalized_MK
#from e45_Sim import statespace,calcOutput
from HW2_library import create_plot
from HW4_library import straight_beam_inertia, drive_train_inertia

#this initial format is just to organize outputs without running everything every time,
#feel free to change

#set to 1 to run all problems
run_problem = 4
debug=True

if run_problem==2 or run_problem==1:
    print(" ")  
    print("Problem 2:")
    #Problem 2  -----------
    
    #import 22MW tower data
    ED = FASTInputFile('ElastoDyn.dat')
    df = FASTInputFile('ElastoDyn_tower.dat').toDataFrame()
    #print(df.keys())

    EI = df["TwFAStif_[Nm^2]"].to_numpy() #stiffness of tower along height
    HtFract = df["HtFract_[-]"].to_numpy() #height fraction along height of tower
    if debug:
        #create_plot(x,y,title,x_axis,y_axis,labels,legend,line):
            #x,y are lists of lists for data sets
            #labels is list of string names for each data set
            #legend,line are a boolean
        create_plot([HtFract],[EI],"Tower Stiffness from ElastoDyn_tower.dat","Height fraction","EI (Nm^2)","EI",False,True)
    
    #gather data to calculate 2 beam modes for tower:
    n_modes = 2
    L_t = ED["TowerHt"]- ED["TowerBsHt"]
    print(f"Tower Height = {L_t:,.5} m" )
    m_prime = df["TMassDen_[kg/m]"].to_numpy() #mass per unit length
   
    #distance along tower height
    z = L_t*HtFract #same number of points as HtFract
    
    freqs, modes = compute_modes(n_modes, z, EI, m_prime)
    freq1 = freqs[0]
    freq2 = freqs[1]
    
    print(f"Beam Mode 1 frequency = {freq1:,.3} Hz " )
    print(f"Beam Mode 2 frequency = {freq2:,.3} Hz " )
    
    #plot mode shapes:
    plt.figure()
    plt.plot(z, modes[0], label="Beam Mode 1")
    plt.plot(z, modes[1], label="Beam Mode 2")
    plt.xlabel("Distance z Along Tower (ground up) (m)")
    plt.ylabel("displacement")
    plt.title("Beam Mode Shapes")
    plt.legend()
    plt.show()
        
    print(" ")   
    #problems 2 and 3 are linked
    print("Problem 3:")
    #Problem 3  -----------
    #calculate and plot tower 22MW tower deflection with tip force
    
    Tmax = 3*(10**6) #Thrust (tip force) in N (3 MN)
    
    #fore-aft deflection
    u, theta, kappa, S, M = deflection(z,EI,Ftip=Tmax)
    
    #plot deflection
    plt.figure()
    plt.plot(z, u, color = "red")
    plt.xlabel("Distance z Along Tower (ground up) (m)")
    plt.ylabel("Deflection (m)")
    plt.ylim(-0.2, 1.5)
    plt.title("Fore-aft Tower Deflection with Tmax = 3MN at tip")
    plt.show()

    umax = u[-1]
    print(f"Max deflection at tower tip = {umax:,.4} m " )
    print("(applied force = 3 MN )")

if run_problem==4 or run_problem==1:
    print(" ") 
    print("Problem 4:")
    #Problem 4  -----------
    
    #part a)
    #see function straight_beam_inertia() in HW4_library.py
    #import 22MW tower data
    ED = FASTInputFile('ElastoDyn.dat')
    df = FASTInputFile('ElastoDyn_tower.dat').toDataFrame()
    L_t = ED["TowerHt"]- ED["TowerBsHt"]
    HtFract = df["HtFract_[-]"].to_numpy()
    z = L_t*HtFract #same number of points as HtFract
    m_prime = df["TMassDen_[kg/m]"].to_numpy() #mass per unit length
    
    #Smart Test Case
    #Assuming a Test tower with constant mass distribution, Height 50m with 
    #m' = 100 kg/m
    
    #analytical test results:
    M_test = 5000 #kg
    S_test = 125000 #kgm
    J_test = 4.1666 * 10**6 #kgm
    
    z_test = np.linspace(0,50,20) #m
    m_test = np.ones(20)*100 #kg/m
    M_t, S_t, J_t = straight_beam_inertia(z_test, m_test)
    
    print("Test Case:")
    print(f" Analytical M = {M_test} kg, calculated M = {M_t} kg")
    print(f" Analytical S = {S_test} kgm, calculated S = {S_t} kgm")
    print(f" Analytical J = {J_test} kgm^2, calculated J = {J_t} kgm^2")
    
    #part b)  Use function on tower data
    M,S,J = straight_beam_inertia(z,m_prime)
    
    print("")
    print("Actual Tower Values: ")
    print(f"Total Mass (M) = {M:,.4} kg ")
    print(f"First Moment of Inertia (S) = {S:,.4} kgm ")
    print(f"Second Moment of Inertia (J) = {J:,.4} kgm^2 ")
    cm = S/M
    print(f"Center of Mass height = {cm:,.4} m ")
    
    #part c)  Drivetrain inertia
    #from ElastoDyn.dat:
    J_gen = 3117597.5395135996 #kg m^2 - generator inertia on HSS
    J_hub = 1.88e+06 # kg.m^2 - hub inertia
    J_hub = 1882498.8832803261
    n = 1.0 #gear ratio of 1.0 from the doc
    
    #dive train:
    J_DT = drive_train_inertia(n, J_hub, J_gen)
    print("")
    print(f"Drive Train Inertia (J_DT)  = {J_DT:,.7} kg m^2 ")

    #part d)
    EI = df["TwFAStif_[Nm^2]"].to_numpy()
    freqs, mode = compute_modes(1, z, EI, m_prime)
    GM , GK_EI , GK_Np, GK_NM = generalized_MK(z, m_prime, EI, mode[0], p_z=None, M_top=1.05e+06, g=9.81)
    print(f"Generalized Mass = {GM:.7} kg and Stiffness = {GK_EI:.7}")


if run_problem==5 or run_problem==1:
    print(" ") 
    print("Problem 5:")
    #Problem 5  -----------

    
    
    
    
    


