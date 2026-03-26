# Import Statements
import os
import matplotlib
import matplotlib.pyplot as plt
#from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd
from scipy.interpolate import interp1d
from libBEM import aero_coeffs,spanloads,intloads,BEMqs


def q4():
    
    #check directory to make sure the right files are reachable
    #print("CWD =", os.getcwd())
    #print("Files here:", os.listdir())
    
    #PART a) -----------------------------------------------------
    print("Problem 4 ----------------")
    print("Part (a): ")
    #this only works when the working directory is the Q4 folder
    #r, twist, chord = pd.read_csv("../data/AeroBlade_low.csv").values.T
    r, twist, chord = pd.read_csv("../data/AeroBlade_High.csv").values.T
    polar = pd.read_csv("../data/ffa-w3-211_Re10M.csv").values
    fPolar = interp1d(polar[:,0], polar[:,1:], axis=0)
    rhub, R = r[0], r[-1]
    r_full = r + rhub
    uin, uit = np.zeros_like(r), np.zeros_like(r)
    #plotting controls
    pitch = 0 #deg
    U0 = 10 #m/s
    B = 3
    lambda_d = 11 #design value
    Omega = lambda_d*U0/R #rad/s
    df = spanloads(pitch, Omega, U0, r_full, uin, uit, chord, twist, fPolar)
    
    print("Plotting Values for: U0 =",U0,"m/s  lambda =",lambda_d,"  Pitch =",pitch,"deg")
    
    #set to true to test calculations against sample data
    overlay_sample_data = True  
    
    alpha = df["alpha_[deg]"].to_numpy()
    ft = df["ft_[N/m]"].to_numpy()
    fn = df["fn_[N/m]"].to_numpy()
    r = df["r_[m]"].to_numpy()
    mz = df["mz_[Nm/m]"].to_numpy()
    
    #Sample test data provided - s marks sample
    r_s = np.array([5, 25, 45, 65, 85, 105, 125, 145, 165, 185])
    alpha_s = np.array([34.48, 21.31, 16.83, 14.84, 13.73, 13.03, 12.55, 12.21, 11.94, 11.73])
    ft_s = np.array([367, 1496, 1945, 2141, 2120, 2074, 2018, 1960, 1834, 100])

    #Plot: loads on left axis and alpha on right axis
    #looked to online sources for assistance with combined plot code
    fig, ax1 = plt.subplots()

    #Calculated values using file
    ax1.plot(r,fn,color = 'green', label="fn(r) [N/m]")
    ax1.plot(r,ft,color = 'blue', label="ft(r) [N/m]")
    ax1.set_xlabel("r [m]")
    ax1.set_ylabel("Spanwise loads [N/m]")
    ax1.grid(True)

    #Overlay sample ft points
    if overlay_sample_data:
        ax1.plot(r_s, ft_s, linestyle="None", marker="X", color = 'blue', label="Sample ft [N/m]")

    #alpha with right axis labelled
    ax2 = ax1.twinx()
    ax2.plot(r, alpha, linestyle="-", color = 'red', label="alpha(r) [deg]")
    ax2.set_ylabel("Alpha [deg]")

    # Overlay sample alpha
    if overlay_sample_data:
        ax2.plot(r_s, alpha_s, linestyle="None", marker="X", color = 'red', label="Sample alpha [deg]")

    # Combine legends from both axes -- AI helped with this chunk
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")
    plt.title("Spanloads Outputs")
    plt.show()
    
    #Part b ----------------------------------------------------------------
    print("Part (b): ")
    #percent of total power per blade section
    #break the blade into sections based on the difference between each r point - dr
    
    #first break blade into n sections of length dr
    n = len(r)-1
    Pi = []
    Pt = 0 #total power
    section = []
    
    for i in range(n):
        #dr of section
        dr = r[i+1]-r[i]
        #midpoint of section for better accuracy
        mid_r = r[i]+0.5*dr
        #linearly extrapolate ft(r) for midpoint of teach section
        ft_mid = (ft[i+1]+ft[i])/2
        dft = ft_mid*dr
        #power
        dP = B*Omega*mid_r*dft
        Pi.append(dP)
        Pt = Pt+dP
        section.append(i)
        
    Pi = (Pi/Pt)*100
    
    plt.figure()
    plt.plot(section, Pi)
    plt.xlabel("Section number")
    plt.ylabel("percentage of total power")
    plt.title("BET Section Pecentage of Total Power Produced")
    plt.show()
    
    #verify calculations, percentages should add up to 100
    tot = sum(Pi)
    print(f"Sum of section power percentages: {tot:.4f}")
    
    #Part c -------------------------------------------------------------
    print("Part (c): ")
    #using intloads to plot Cp(lambda)
    #from range of 6 to 12
    lambda_d = np.linspace(1,12,100)
    Omega = lambda_d*U0/R
    
    #for pitch = 0 degrees
    pitch = 0
    CP0 = []
    for i in range(len(Omega)):
        df = spanloads(pitch, Omega[i], U0, r_full, uin, uit, chord, twist, fPolar)
        ft = df["ft_[N/m]"].to_numpy()
        fn = df["fn_[N/m]"].to_numpy()
        r = df["r_[m]"].to_numpy()
        mz = df["mz_[Nm/m]"].to_numpy()
        loads = intloads(r, fn, ft, mz, R, U0, Omega[i], rho=1.225, B=3)
        Cp = loads["CP"]
        CP0.append(Cp)
    
    #for pitch = 20 degrees
    pitch = 20
    CP20 = []
    for i in range(len(Omega)):
        df = spanloads(pitch, Omega[i], U0, r_full, uin, uit, chord, twist, fPolar)
        ft = df["ft_[N/m]"].to_numpy()
        fn = df["fn_[N/m]"].to_numpy()
        r = df["r_[m]"].to_numpy()
        mz = df["mz_[Nm/m]"].to_numpy()
        loads = intloads(r, fn, ft, mz, R, U0, Omega[i], rho=1.225, B=3)
        Cp = loads["CP"]
        CP20.append(Cp)
        
    #for pitch = 40 degrees
    pitch = 40
    CP40 = []
    for i in range(len(Omega)):
        df = spanloads(pitch, Omega[i], U0, r_full, uin, uit, chord, twist, fPolar)
        ft = df["ft_[N/m]"].to_numpy()
        fn = df["fn_[N/m]"].to_numpy()
        r = df["r_[m]"].to_numpy()
        mz = df["mz_[Nm/m]"].to_numpy()
        loads = intloads(r, fn, ft, mz, R, U0, Omega[i], rho=1.225, B=3)
        Cp = loads["CP"]
        CP40.append(Cp)
        
    plt.figure()
    plt.plot(lambda_d, CP0, label="Pitch = 0 degrees")
    plt.plot(lambda_d, CP20, label="Pitch = 20 degrees")
    plt.plot(lambda_d, CP40, label="Pitch = 40 degrees")
    plt.xlabel("lambda (TSR)")
    plt.ylabel("Cp")
    plt.title("Cp as a function of Tip Speed Ratio")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # Code here is not executed in main
    q4()
    # Code

