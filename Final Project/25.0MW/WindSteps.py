import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import __editable___openfast_toolbox_3_5_1_finder

# Wind Step Input
from openfast_toolbox.case_generation.case_gen import createStepWind

createStepWind('WindStepPowerCurve3.wnd', WSstep=1, WSmin=2, WSmax=15, tsstep=100, dt=0.5)

# trying to plot
# plt.plot(Cd, Cl, marker='o', color='red')
# plt.plot(cds, cls, marker='o')
# plt.xlabel('Coefficient of Drag (Cd)')
# plt.ylabel('Coefficient of Lift (Cl)')
# plt.title('Total Pressure Force Function')
# plt.show()

# plt.plot(cds, cls, marker='o')
