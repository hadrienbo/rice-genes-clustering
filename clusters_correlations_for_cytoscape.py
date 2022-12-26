import math
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm        # module palettes de couleurs
import pandas as pd
import seaborn as sns

fichier="/content/drive/MyDrive/PFE/meanPoints_CLUSTERING1_R2_pow1.csv"
clusters_data2 = pd.read_csv(fichier,sep=",")
clusters_data2 = clusters_data2.rename(columns = {'Unnamed: 0' : 'Name'})
clusters_data2=clusters_data2.set_index('Name')
clusters_data2
fig=plt.figure(figsize=(22.5,15))
from  matplotlib.colors import LinearSegmentedColormap
c = ["darkred","red","lightcoral","white", "palegreen","green","darkgreen"]
v = [0,.15,.4,.5,0.6,.9,1.]
l = list(zip(v,c))

cmap=LinearSegmentedColormap.from_list('rg',l, N=256)
sns.heatmap(clusters_data.transpose().corr(), annot=True,cmap=cmap)
plt.show()

fig=plt.figure(figsize=(22.5,15))
cmap2=LinearSegmentedColormap.from_list('rg',l, N=256)
sns.heatmap(clusters_data2.transpose().corr(), annot=True,cmap=cmap2)
plt.show()

fichier="/content/drive/MyDrive/PFE/new_SFs.csv"
epigenetic_data = pd.read_csv(fichier,sep=",")
epigenetic_data=epigenetic_data.set_index('Name')
epigenetic_data.columns=['X0h.brep1',	'X0h.brep2', 'X0h.brep3',	'X1h.brep1',	'X1h.brep2',	'X1h.brep3',	'X3h.brep1',	'X3h.brep2',	'X3h.brep3',	'X6h.brep1',	'X6h.brep2',	'X6h.brep3',	'X12h.brep1',	'X12h.brep2',	'X12h.brep3',	'X24h.brep1',	'X24h.brep2',	'X24h.brep3',	'X36h.brep1',	'X36h.brep2',	'X36h.brep3',	'X48h.brep1',	'X48h.brep2',	'X48h.brep3']

df=pd.concat([clusters_data2,epigenetic_data],axis=0)
CORRELATIONS=df.transpose().corr().iloc[:20,20:]
CORRELATIONS.to_csv('/content/drive/MyDrive/PFE/NEW_SFandCLUSTconnexions.csv')
fig=plt.figure(figsize=(22.5,7.5))
sns.heatmap(CORRELATIONS, annot=False)
plt.show()

flat_corr=pd.DataFrame(CORRELATIONS.to_numpy().flatten())
flat_corr.nlargest(38,columns=0)
CORRELATIONS.where(abs(CORRELATIONS)>0.915899)
input=[]
output=[]
corr=[]
for row in CORRELATIONS.index:
  for col in CORRELATIONS.columns:
    if CORRELATIONS.loc[row,col]>0.915899 :
      input.append(col)
      output.append(row)
      corr.append(CORRELATIONS.loc[row,col])
connexions=pd.DataFrame(np.array([input,output,corr]).T)
connexions.columns=['Splic','Cluster','Correlation']
fichier="/content/drive/MyDrive/PFE/connexions_1%_Epig+R2_69_M0.2.csv"
connexions2 = pd.read_csv(fichier,sep=",")
connexions2 = connexions2.rename(columns = {'Unnamed: 0' : 'Name'})
connexions2=connexions2.set_index('Name')
