import os
import matplotlib
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
import numpy as np
import math
import pandas as pd

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

U_avg = np.mean(avg_speed_mph)
print(f"Yearly Average Wind Speed = {U_avg:,.4} mph " )