# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 15:26:13 2026

@author: Owner
"""

import os
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np

def create_plot(x,y,title,x_axis,y_axis,labels,legend,line):
    #my own sourcecode
    #x,y are lists of lists for data sets
    #labels is list of string names for each data set
    #legend,line are a boolean
    
    fig = plt.figure() #initiate plot
    
    #Loop through data sets to plot
    for i in range(len(x)):
        if line==True:
            plt.plot(x[i],y[i],label=labels[i])
        else:
            plt.plot(x[i],y[i],label=labels[i],linestyle="None", marker="o")
        
    if legend == True:
        plt.legend()
        
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.title(title)
    plt.show()

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

def freq_peaks(fa,S_ka,n):
    #returns the n number of frequency peaks in a spectrum
    
    fa = np.asarray(fa)
    S_ka = np.asarray(S_ka)

    if fa.size == 0 or S_ka.size == 0:
        raise ValueError("Empty input arrays.")
    if fa.size != S_ka.size:
        raise ValueError(f"Length mismatch: len(fa)={fa.size}, len(S_ka)={S_ka.size}")

    # --- n == 1: global max ---
    if n == 1:
        idx = int(np.nanargmax(S_ka))   # ignores NaNs; use argmax if no NaNs
        return float(S_ka[idx]), float(fa[idx])

    # --- find local maxima (strict) ---
    # peak if S[i-1] < S[i] > S[i+1]
    peaks = np.where((S_ka[1:-1] > S_ka[:-2]) & (S_ka[1:-1] > S_ka[2:]))[0] + 1

    if peaks.size == 0:
        return [], []

    # take top-n peaks by PSD magnitude
    order = np.argsort(S_ka[peaks])[::-1]
    top = peaks[order[:n]]

    return list(S_ka[top].astype(float)), list(fa[top].astype(float))
    
    """
    #first convert array to list
    S_k = []
    f = []
    for i in range(len(S_ka)):
        S_k[i] = S_ka[i]
        f[i] = fa[i]
    S_list = []
    f_list = []
    if n==1:
        return max(S_k),f[S_k.index(max(S_k))]
    else:
        for i in range(len(f)):
            if i>1:
                if S_k[i]>S_k[i-1] and S_k[i]>S_k[i:i+2]:
                    S_list.append(S_k[i])
                    f_list.append(f[i])
        sl = []
        fl = []
        for i in range(n):
            sl.append(max(S_list))
            fl.append(f_list[S_list.index(max(S_list))])
            del f_list[S_list.index(max(S_list))]
            del S_list[S_list.index(max(S_list))]
    return sl,fl
    """

def main():
    
    """
    x1 = [[1,2,3,4,5],[1,2,4,9,16,25]]
    y1 = [[1,4,6,8,10],[1,2,3,4,5,6]]
    x2 = [[1,10,15,23,40]]
    y2 = [[5,3,8,54,2]]
    L1 = ["set 1","set 2"]
    L2 = ["AHH"]
    
    create_plot(x1,y1,'Data set 1','X','Y',L1,True)
    create_plot(x2,y2,'Data set 2','X','Y',L2,False)
    """
    
if __name__ == "__main__":
    main()