import math
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm        # module palettes de couleurs
import pandas as pd
import seaborn as sn

fichier="/content/drive/MyDrive/PFE/TPM_trans.csv"
df3 = pd.read_csv(fichier,sep=",",index_col=0)
df3=df3.loc[:,['0h.brep1', '0h.brep2', '0h.brep3','1h.brep1', '1h.brep2', '1h.brep3',
       '3h.brep1', '3h.brep2', '3h.brep3','6h.brep1','6h.brep2', '6h.brep3', '12h.brep1', '12h.brep2', '12h.brep3', '24h.brep1', '24h.brep2', '24h.brep3', '36h.brep1', '36h.brep2',
       '36h.brep3', '48h.brep1', '48h.brep2', '48h.brep3']]
files=[df3]
for f in files:
  for i in range(f.shape[0]):
    print(i/f.shape[0])
    m=np.nanmean(np.array(f.iloc[i,:]))
    std=np.std(np.array(f.iloc[i,:]))
    for j in range(f.shape[1]):
      f.iloc[i,j]=(f.iloc[i,j]-m)/std
      
