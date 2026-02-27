from libwind import kaimal_f , NTM_TI
from libsignal import psd , autoCorrCoeff , bin_DF
import pandas as pd
import glob

# Find all .txt files in the current directory
files = glob.glob("*.parquet")
print(files)
lv = 0
df = []
for f in files:
    df.append(pd.read_parquet(f))
    lv += 1
print(df[0])


