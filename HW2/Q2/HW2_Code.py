#Nate Fay -- HW2 -- Problems 1 and 3
import os
import matplotlib
import matplotlib.pyplot as plt
#from openfast_toolbox.io import FASTInputFile
import numpy as np
#matplotlib.use("TkAgg")  # or "QtAgg"

from libwind import kaimal_f , NTM_TI
from libsignal import psd , autoCorrCoeff , bin_DF
from HW2_library import create_plot, loglog_slope, freq_peaks
import pandas as pd
from scipy.stats import weibull_min
import math

def main():
   
    #note: this code runs the calculations and plots for problem 1 and 3 of HW2
    #use the top variable p to toggle between problem outputs
    #sections that are commented out are code that I attempted but had errors or bugs,
    #that I replaced with improved versions that accomplish the same goal with the help of AI
    #the accompanying files to have in the same directory are:
        #HW2_library.py, libwind and libsignal
   
    
    p = 3 #select problem number to display outputs/plots
    #1 or 3 only for this code file
   
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
    
    #problem 2
    
    
    #problem 3
    if p==3:
        print("part A")
        filepath = 'C:/Users/Owner/Documents/MIE673/HW2/HW2_Wind/HW2_M2Data_sel2.parquet'
        df = pd.read_parquet(filepath) 
        WS80 = df['WS_80_[m/s]']
        WD80 = df['WD_80_[deg]']
        sig80 = df['sigma_80_[m/s]']
        
        #auto weibull fit method
        wind = np.asarray(WS80)
        # Fit Weibull 
        k, loc, lam = weibull_min.fit(wind, floc=0)
        x = np.linspace(0, wind.max(), 500)
        pdf = weibull_min.pdf(x, k, loc=0, scale=lam)

        #weibull fit based on Wind Energy Explained formula
        #weibull parameters
        hist, bin_edges = np.histogram(wind, bins = 40)
        dataAverage = np.average(WS80)                      # find average
        dataStdev = np.std(WS80)
        k2 = (dataStdev/dataAverage)**(-1.086)
        c = dataAverage/(math.gamma(1+(1/k)))
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        U = np.zeros(40) # initialize wind speed
        hist_2 = np.zeros(40) # initialize normalized histogram
        pdf2 = np.zeros(40) # initialize pdf
        for i in range(40):
            U = bin_centers[i] # load wind speeds into U
            pdf2[i] = (k2/c)*(U/c)**(k2-1)*math.exp(-(U/c)**k2) # create the pdf

        
        plt.figure()
        plt.hist(wind, bins=40, density=True, label="Data (PDF)")
        plt.plot(x, pdf, label=f"Weibull fit (k={k:.2f}, A={lam:.2f})")
        U = bin_centers
        #plt.plot(U, pdf2, linewidth=2, color = 'red', label='WEE 3E Weibull PDF')
        plt.xlabel("Wind speed [m/s]")
        plt.ylabel("Probability density [1/(m/s)]")
        plt.title("Wind speed at 80m PDF with Weibull fit")
        plt.legend()
        plt.grid(True)
        plt.show()
        
        print("WEE3E Weibull c: ", c, " k: ", k2)
        
        #calculate hours per year spent above 25 m/s
        #sampling rate = 10 min
        #by hour, fs = 6 Hz
        #for 2 years)
        #using weibull parameters
        v0 = 25.0
        p_exceed = np.exp(-(v0/lam)**k)
        hours_per_year = 8760 * p_exceed
        print("P(V>25) =", p_exceed)
        print("Expected hours/year >25 m/s =", hours_per_year)
        
        #from data
        v = np.asarray(WS80)
        p_exceed_emp = np.mean(v > v0)
        #p_exceed_emp = len(v[v>v0])/len(v)
        hours_per_year_emp = 8760 * p_exceed_emp
        print("Empirical P(V>25) =", p_exceed_emp)
        print("Empirical hours/year >25 m/s =", hours_per_year_emp)
        
        #part b
        print("part B")
        U = WS80
        U_turbulence = U - np.mean(U) #mean removed to get turbulence
        time = []
        for i in range(len(U_turbulence)):
            time.append(i*600) #convert to seconds (10 min samples)
        #f, S_xf = psd(time,U_turbulence,'welch',2**10)
        f, S_xf = psd(time,U_turbulence,None,None)
        #create_plot([f],[S_xf],"Wind Data 80m","Frequency (Hz)","PSD",["data1"],False,True)
       
        '''#original code below: had to be improved by AI to remove error
        n = 20 #number of peaks
        Sx, fx = freq_peaks(f, S_xf, n )
        for i in range(n):
            print("Peak Frequency: ", fx[i], "Hz")
        plt.figure()
        plt.plot(f,S_xf,label="psd of data")
        for i in range(n):
            plt.plot(fx[i],Sx[i],color='red',marker='o', markersize=1)
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("PSD (m/s)^2/Hz")
        plt.title("Wind Data 80m")
        plt.show
        '''
        n = 20 #number of frequency peaks to identify
        #identify and print peak values
        Sx, fx = freq_peaks(f, S_xf, n)

        n_found = min(len(fx), len(Sx))
        print("peaks requested:", n, "peaks found:", n_found)
        
        for i in range(n_found):
            print("Peak Frequency:", fx[i], "Hz")

        plt.figure()
        plt.plot(f, S_xf, label="psd of data")

        # plot peaks (bigger marker so you can see them)
        plt.scatter(fx[:n_found], Sx[:n_found], s=5, color="red", zorder=3)

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("PSD (m/s)^2/Hz")
        plt.title("Wind Data 80m")
        plt.show(block=True)
        
        #part c
        print("part C")
        u = np.asarray(WS80)  #wind speed time series [m/s]
        bin_width = 1.0
        bins = np.arange(0, np.nanmax(u) + bin_width, bin_width)
        sigma = []
        u_mid = []

        for i in range(len(bins) - 1):
            lo, hi = bins[i], bins[i+1]
            mask = (u >= lo) & (u < hi)

            if np.sum(mask) >= 5:  # require a few points
                sigma.append(np.std(u[mask], ddof=0))
                u_mid.append(0.5 * (lo + hi))
            else:
                sigma.append(np.nan)
                u_mid.append(0.5 * (lo + hi))

        sigma = np.array(sigma)
        u_mid = np.array(u_mid)

        plt.figure()
        plt.plot(u_mid, sigma, marker="o", linestyle="None")
        plt.xlabel("Wind speed U (bin midpoint) [m/s]")
        plt.ylabel("Standard deviation σ(U) [m/s]")
        plt.title("σ as a function of wind speed with 1 m/s bins")
        plt.grid(True)
        plt.show()
        
        #compute TI of each bin
        '''
        df_TI = []
        dfb_mean = bin_DF(df , colBin='WS_80_[m/s]', xBins=range(0, 35, 1),stats='mean')
        dfb_std = bin_DF(df , colBin='WS_80_[m/s]', xBins=range(0, 35, 1),stats='std')
        
        #now create data set of TI using the formula
        for i in range(len(dfb_mean)):
            df_TI.append(dfb_std[i]/dfb_mean[i])
        
        plt.figure()
        plt.plot(dfb_mean, dfb_std, marker="o", label='STDV scatter',linestyle="None")
        plt.plot(dfb_mean, df_TI, marker="o", label='TI')
        plt.xlabel("Wind speed U (bins) [m/s]")
        plt.ylabel("Standard deviation σ(U) [m/s] or Turbulence Intensity")
        plt.title("Turbulence Intensity vs Wind Speed 80m")
        plt.grid(True)
        plt.show()
        '''
        
        xBins = np.arange(0, 36, 1)

        dfb_mean = bin_DF(df.copy(), colBin='WS_80_[m/s]', xBins=xBins, stats='mean')
        dfb_std  = bin_DF(df.copy(), colBin='WS_80_[m/s]', xBins=xBins, stats='std')

        # x: mean wind speed in each bin
        U_bin = dfb_mean['binMean'].to_numpy()

        # y: std in each bin (column has same name as colBin)
        sigma_bin = dfb_std['WS_80_[m/s]'].to_numpy()

        TI_bin = sigma_bin / U_bin
        TI_bin[U_bin <= 0] = np.nan  # avoid divide-by-zero

        #IEC models
        modelA_TI = NTM_TI(U_bin,'A')
        modelB_TI = NTM_TI(U_bin,'B')
        modelC_TI = NTM_TI(U_bin,'C')

        plt.figure()
        plt.plot(U_bin, sigma_bin, marker="o", linestyle="None", label="σ(U) binned")
        plt.plot(U_bin, TI_bin, marker="o", linestyle="-", label="TI(U) binned")
        plt.plot(U_bin, modelA_TI, linestyle="-", label="IEC Class A model")
        plt.plot(U_bin, modelB_TI, linestyle="-", label="IEC Class B model")
        plt.plot(U_bin, modelC_TI, linestyle="-", label="IEC Class C model")
        plt.xlabel("Wind speed U [m/s]")
        plt.ylabel("σ(U) [m/s]  and  TI (unitless)")
        plt.title("Binned σ(U) scatter and TI vs Wind Speed (80 m)")
        plt.grid(True)
        #plt.ylim(0, 1)
        plt.legend()
        plt.show()
        print("nothing to print, see plots")
        
        #part d - power law
        print("part D")
        WS80 = df['WS_80_[m/s]']
        WS50 = df['WS_50_[m/s]']
        WS20 = df['WS_20_[m/s]']
        WS10 = df['WS_10_[m/s]']
        WS5 = df['WS_5_[m/s]']
        WS2 = df['WS_2_[m/s]']
        
        #calculate mean wind speed at each height
        U80 = WS80.mean()
        U50 = WS50.mean()
        U20 = WS20.mean()
        U10 = WS10.mean()
        U5 = WS5.mean()
        U2 = WS2.mean()
        #mean_U = [U2,U5,U10,U20,U50,U80]
        #heights = [2,5,10,20,50,80]
        #convert to array to use better functions
        mean_U = np.array([U2, U5, U10, U20, U50, U80], dtype=float)
        heights = np.array([2, 5, 10, 20, 50, 80], dtype=float)

        #Use only positive values
        mask = (mean_U > 0) & (heights > 0)
        x = np.log(heights[mask])
        y = np.log(mean_U[mask])
        alpha, intercept = np.polyfit(x, y, 1)  #y = alpha*x + intercept
        print("model fit alpha =", alpha)
        #now create model to plot against the data to compare
        z_ref = 10.0 #arbitrary reference point for the model
        U_ref = U10
        h = np.linspace(0.1,85,500)
        U_model = U_ref * (h / z_ref)**alpha
        plt.figure()
        plt.plot(mean_U, heights, "o", label="Measured means")      # points
        plt.plot(U_model, h, "-", label=f"Power law model fit (α={alpha:.3f})")  # model curve
        plt.xlabel("Wind speed (m/s)")
        plt.ylabel("Height (m)")
        plt.title("Mean windspeed vs power-law model")
        plt.grid(True)
        plt.legend()
        plt.show()
           
if __name__ == "__main__":
    main()



