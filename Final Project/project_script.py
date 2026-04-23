import os
import matplotlib
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd
import project_library.py as lib

#wind
import numpy as np

months = [
    "January","February","March","April","May","June",
    "July","August","September","October","November","December","Annual"
]

# Normal Monthly Average Speed [MPH]
avg_speed_mph = np.array([
    45.6, 44.6, 39.8, 35.6, 29.6, 26.8,
    25.5, 24.5, 27.6, 35.5, 39.4, 44.0, 34.9
], dtype=float)
print(f"Yearly Average Wind Speed = {U_avg:,.4} mph " )

U_avg = np.mean(avg_speed_mph)
C_P_rated_0 = 0.3
P = 25 * 10**6
rho = 1.225 # wrong but it's fine
R_0 = lib.findR(P, rho, U, C_P)
print(R_0)
