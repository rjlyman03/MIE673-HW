# Setup path
import sys
sys.path.append("../mie673-data/sample_codes")

import os
from glob import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTInputFile
from e01_PlotBlade import blade_chord
# 
wd = os.getcwd()
data_root ='../mie673-data'
folders = glob("**MW", root_dir = data_root)
N=len(folders)
cm = mpl.colormaps['viridis']
colors = cm(np.linspace(0, 1, N))
fig, axs = plt.subplots(N, sharex=True, sharey=True, layout="constrained")
fig.suptitle('Chord Length vs Span for different MW Wind turbines', fontweight='bold')
for n in range(N): 
    case = folders[n]
    os.chdir(data_root + '/' + case)
    r, chord = blade_chord('AeroDyn_blade.dat')

    # --- Plot
    axs[n].plot(r, chord, color=colors[n])
    axs[n].set_title(case, pad=12) 
    axs[n].grid(True)
    axs[n].set_yticks([0, 2, 4, 6, 8])
    #if n != (N-1):
    #    axs[n].set_xticks([])
    os.chdir(wd)

fig.supxlabel('Span [m]', fontweight='bold')
fig.supylabel('Chord [m]', fontweight='bold')
plt.show()
fig.savefig('new.png')
