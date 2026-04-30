import pandas as pd
import numpy as np
import aerodesign as ad
import libBEM as bem
import project_library as lib

df = lib.readPolarDatFile("./22.0MW/Airfoils/IEA-22-280-RWT_AeroDyn15_Polar_56.dat")
cd = np.asarray([df["Cd"], df["Cl"]])
print(cd)