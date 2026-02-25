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