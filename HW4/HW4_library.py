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

def straight_beam_inertia(z,m):
    #z is an array of distances (m)
    #m is an array of linear mass densitys corresponding to z (kg/m)
    M = np.trapezoid(m, z)
    S = np.trapezoid(m*z, z)
    J = np.trapezoid(m*(z**2))
    return M , S , J
