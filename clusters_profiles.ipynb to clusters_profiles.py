import matplotlib
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 28}

matplotlib.rc('font', **font)
import math
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm        # module palettes de couleurs
import pandas as pd
import seaborn as sns

plt.rcParams['font.size'] = 18
plt.rcParams['figure.figsize'] = (18, 6)
#np.set_printoptions(edgeitems=10)
np.set_printoptions(linewidth = 220)
np.set_printoptions(precision=4)
#np.set_printoptions(precision=3,formatter={'float': '{:9.3f}'.format})
pd.set_option('precision', 5)
pd.set_option("display.max_columns",20)
pd.set_option('display.max_rows', 999)
#pd.set_option('max_colwidth', 6)
linestyles = [(0, ()), # solid 
              (0, (5, 10)),(0, (5, 5)),(0, (5, 1)), # dashed (loosely/normal/densely)
              (0, (3, 10, 1, 10)),(0, (3, 5, 1, 5)),(0, (3, 1, 1, 1)), # dotted  (loosely/normal/densely)
              (0, (3, 10, 1, 10, 1, 10)),(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1, 1, 1)), # dashdotted(loosely/normal/densely)
              (0, (1, 10)),(0, (1, 5)),(0, (1, 1))]
couleurs = cm.Dark2.colors
fichier="/content/drive/MyDrive/PFE/T_DTU_Zscore.csv"
data = pd.read_csv(fichier,sep=",")
fichier="/content/drive/MyDrive/PFE/CLUSTERING_hclustKmeans.csv"
#fichier="/content/drive/MyDrive/PFE/Clustering1_normalized_data.csv"
clustering = pd.read_csv(fichier,sep=",")
clustering = clustering.rename(columns = {'Unnamed: 0' : 'Name'})
clustering=clustering.set_index('Name')
clusters_mean_points=[]
nb_clusters=max(clustering.iloc[:,2])
N=len(clustering)
print("nombre de clusters =",nb_clusters+1)
#fig, axs = plt.subplots(nb_clusters)
#fig.suptitle('Clusters profiles')
plt.figure()
X = data.to_numpy().T
X=X[1:,:]
nomDesVariables=list(data.iloc[:,0])
nomDesIndividus=list(data.columns[1:])
fig, ax = plt.subplots(nrows=int((nb_clusters+1)/2), ncols=1, figsize=(50,150))                          
cluster=10
for row in ax:
  cluster_indexes=[]
  for i in range(N):
    if clustering.iloc[i,2]==cluster:
      cluster_indexes.append(i)
  X0=X[cluster_indexes,:]
  moyennes = X0.mean(axis=0)
  clusters_mean_points.append(moyennes)
  #label_variables="TPM"
  label_variables="Z-score"
  for i, (ligne,label) in enumerate(zip(X0, nomDesIndividus)):
    row.plot(ligne, label=label,
             color = couleurs[i%len(couleurs)],
             linestyle=linestyles[(i//len(couleurs))%len(linestyles)])
  row.plot(moyennes,'k-',label="moyenne",linewidth=3)
  row.title.set_size(70)
