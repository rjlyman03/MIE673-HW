# Import Statements
import os
import matplotlib
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd
from libbeam import deflection,compute_modes,generalized_MK
from e45_Sim import statespace,calcOutput
from HW2_library import create_plot

#this initial format is just to organize outputs without running everything every time,
#feel free to change

#set to 1 to run all problems
run_problem = 1
debug=True

if run_problem==2 or run_problem==1:
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
    
    #plot mode shapes:
    plt.figure()
    plt.plot(z, modes[0], label="Beam Mode 1")
    plt.plot(z, modes[1], label="Beam Mode 2")
    plt.xlabel("distance z along tower (ground up)")
    plt.ylabel("displacement")
    plt.title("Beam Mode Shapes")
    plt.legend()
    plt.show()
        
        
        
        
    print(" ")   
if run_problem==3 or run_problem==1:
    print("Problem 3:")
    #Problem 3  -----------
    



if run_problem==4 or run_problem==1:
    print("Problem 4:")
    #Problem 4  -----------
    
    

if run_problem==5 or run_problem==1:
    print("Problem 5:")
    #Problem 5  -----------

    
    
    
    
    


