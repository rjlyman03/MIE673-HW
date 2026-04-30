import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from openfast_toolbox.io import FASTInputFile

bld = FASTInputFile('AeroDyn_blade.dat')
# print(bld)
print(bld.keys())

nr = len(bld['BldAeroNodes'])
print('Number of radial positions:', nr)

# bld['NumBlNds'] =  nr


R = 120
r_new = np.linspace(0, R, nr)

# Call planform function
c_new     = r_new * 1
twist_new = r_new * 0 + 10



bld['BldAeroNodes'] = np.zeros((nr, 7))

#  0:   BlSpn    1: BlCrvAC ->0   2:BlSwpAC->0 3: BlCrvAng->0    4->BlTwist  5:BlChord          BlAFID
bld['BldAeroNodes'][:,0] =  r_new   # r
bld['BldAeroNodes'][:,1] =  r_new*0 #
bld['BldAeroNodes'][:,2] =  r_new*0 #
bld['BldAeroNodes'][:,3] =  r_new*0 #
bld['BldAeroNodes'][:,4] = twist_new # twsit
bld['BldAeroNodes'][:,5] = c_new # chord
bld['BldAeroNodes'][:,6] =  1 # ID


bld.write('AeroDyn_blade_out.dat')



if __name__ == '__main__':
    pass
