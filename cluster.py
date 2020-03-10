#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np
import umap
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture as GMM



def generate_data(n_spots : int,
                  n_genes : int,
                  n_types : int,
                  )->np.ndarray:

    probs = np.random.uniform(0,1, size = n_genes).reshape(1,-1)
    data = np.zeros((1,n_genes))
    props = np.random.dirichlet(1 * np.ones(n_types))
    nums = np.ceil((n_spots * props)).astype(int)
    pos = [0]

    for z in range(n_types):
        rates  = np.random.uniform(2,
                                   15,
                                   size = n_genes).reshape(1,-1)

        datum = np.random.negative_binomial(rates,
                                            probs,
                                            size = (nums[z],n_genes))
        pos.append(datum.shape[0])
        data = np.vstack((data,datum))

    print(data)

    data = data[1::,:]
    pos = np.array(pos)
    pos = np.cumsum(pos)
    return (data,pos)


np.random.seed(69)
pca = PCA(n_components = 20)
dimred = umap.UMAP()

n_types = 10

x,idx = generate_data(5000,500,n_types)
x = x / x.sum(axis = 1).reshape(-1,1)
z = pca.fit_transform(x)
y = dimred.fit_transform(z)

gmm = GMM(n_components = n_types)
clabs = gmm.fit_predict(X=y,)
vals = np.zeros(x.shape[0])

for pos in range(idx.shape[0]-1):
    vals[idx[pos]:idx[pos+1]] = pos

fig, ax = plt.subplots(1,1)
ax.scatter(y[:,0],
           y[:,1],
           c = vals,
           edgecolor = 'black',
           alpha = 0.7,
           cmap = plt.cm.jet)

ax.set_aspect('equal')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xticks([])
ax.set_yticks([])

ax.set_xlabel('UMAP-1',fontsize = 25)
ax.set_ylabel('UMAP-2',fontsize = 25)

for pos in ['top','right']:
    ax.spines[pos].set_visible(False)
for pos in ['left','bottom']:
    ax.spines[pos].set_lw(3)

plt.show()
