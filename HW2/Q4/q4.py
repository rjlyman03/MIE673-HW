""" Homework 2 - Question 4 """
import matplotlib.pyplot as plt
import numpy as np
from pydatview.tools.spectral import fft_wrap
import pydatview.io as weio

# --- Generated TurbSim Data
# turbsim TurbSimInput1.dat # hub height of 180 m and 10 by 10 YZ grid
# turbsim TurbSimInput2.dat # hub height of 90 m and 21 by 21 YZ grid

# --- TurbSim Output Files
filenames = []
filenames += ['/home/rjlyman03/Documents/2_Schoolwork/Graduate/MIE673/HW/HW2/TurbSimInput1.bts']
filenames += ['/home/rjlyman03/Documents/2_Schoolwork/Graduate/MIE673/HW/HW2/TurbSimInput2.bts']

# --- Data for different actions

# --- Open and convert files to DataFrames
dfs = {}
for iFile, filename in enumerate(filenames):
    dfs_or_df = weio.read(filename).toDataFrame()
    if isinstance(dfs_or_df, dict):
        for k,df in dfs_or_df.items():
            dfs[k+f'{iFile}'] = df
    else:
        dfs[f'tab{iFile}'] = dfs_or_df

# --- Plot preprocessing
def preproXY(x, y):
    x, y, Info = fft_wrap(x, y, output_type='PSD', averaging='Welch', averaging_window='Hamming', detrend=False, nExp=11, nPerDecade=11)
    return x, y

# --- Simple Plots to show PSD for kaimal process
idx = 0
fig, ax = zip(*[plt.subplots() for _ in range(2)])
titles = ['TurbSim Box 1', 'TurbSim Box 2']
for tabName in ['ZMidLine0', 'ZMidLine1']:
    fig[idx].subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    df =  dfs[tabName]
    ax[idx].plot(*preproXY(df['t_[s]'], df['u_[m/s]']), ls='-', lw=1.5, label='u [m/s]')
    ax[idx].set_title(titles[idx])
    ax[idx].set_xlabel('Frequency [Hz]')
    ax[idx].set_ylabel('PSD(u) [(m/s)^2/s]')
    ax[idx].grid()
    ax[idx].set_xscale('log', nonpositive='clip')
    ax[idx].set_yscale('log', nonpositive='clip')
    idx += 1
plt.show()
