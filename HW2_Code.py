
import os
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np

from libwind import kaimal_f , NTM_TI
from libsignal import psd , autoCorrCoeff , bin_DF
from HW2_library import create_plot
import pandas as pd

#I created this function with the help of AI
def loglog_slope(f, S, fmin, fmax):
    f = np.asarray(f)
    S = np.asarray(S)

    # select region and keep positive values only
    m = (f >= fmin) & (f <= fmax) & (f > 0) & (S > 0)

    x = np.log10(f[m])
    y = np.log10(S[m])

    # y = a*x + b  -> a is the log-log slope
    a, b = np.polyfit(x, y, 1)
    return a, b

def main():
   
    p = 1 #select problem number to display outputs/plots
   
    #change filepath below
    #H2.1 section
    filepath = 'C:/Users/Owner/Documents/MIE673/HW2/HW2_Wind/HW2_Kaimal_10mps_I15_L300_3600s.parquet'
    df = pd.read_parquet(filepath) 
    time = df["Time_[s]"]
    U = df["WindSpeed_[m/s]"]
    if p==1:
        create_plot([time],[U],"Wind Data","Time (s)","U (m/s)",["data1"],False,True)
        #create_plot(x,y,title,x_axis,y_axis,labels,legend,line)
   
    U_turbulence = U - np.mean(U) #mean removed to get turbulence
    f, S_xf = psd(time,U_turbulence,'welch',2**10)
   
    #Kaimal distribution parameters
    #3600s at 100 Hz
    U0 = 10 # reference velocity (m/s)
    L = 300 #length scale (m)
    I = 0.15 #Turbulence Intensity (15%)
    stdv = I*U0 #standard deviation (m/s)
    #kaimal_f(f, U0, sigma, L)
    S_k = kaimal_f(f,U0,stdv,L)
    
    if p==1:
        plt.figure()
        plt.loglog(f, S_xf,label="psd of data")
        plt.loglog(f, S_k,label="Kaimal distribution")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("PSD (m/s)^2/Hz")
        plt.title("Turbulence PSD")
        plt.grid(True, which="both")
        plt.legend()
        plt.show()
    
    #compare variance of time series
    #variance of time series
    variance = np.var(U_turbulence)
    #Integrate PSD over frequency 
    #first recalculate S_xf without the scaaling factor used in psd() (took me a long time to figure out the error)
    f, S_xf = psd(time,U_turbulence)
    psd_area = np.trapezoid(S_xf, f)   #trapezoidal integral
    if p==1:
        print("Variance (time series):", variance)
        print("Integral of PSD:", psd_area)
        print("Ratio PSD variance/time series variance:", psd_area / variance)
        #calculate slope between the range of 0.1Hz to 40 Hz
        slope, intercept = loglog_slope(f, S_xf, fmin=0.1, fmax=40.0)
        print("log-log slope =", slope)
        print("-5/3 = ", (-5/3))
    
    #part b - autocorrelation
    
    #takes the AC for 10000 time delays with a step of 1 for the first 10 minutes
    #autoCorrCoeff(x, nMax=None, dt=1, method='corrcoef'):
    rho,tau = autoCorrCoeff(U[:60000], 50000, 1, 'corrcoef')
    create_plot([tau],[rho],"Autocorrelation of Wind Data","Time step (s)","Autocorrelation Coefficient",["data1"],False,True)
    
if __name__ == "__main__":
    main()



