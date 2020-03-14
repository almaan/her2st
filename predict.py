#!/usr/bin/env python3


import pandas as pd
import numpy as np


pth = "/home/alma/w-projects/spatential/data/breast-cancer/curated/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.tsv.gz"

c_pth = "/tmp/tls/tls-associated.tsv"

#ncnt = pd.read_csv(pth,sep = '\t',header = 0, index_col = 0)

# cnt = pd.DataFrame(cnt.values / cnt.values.sum(axis = 1,keepdims = True),
#                    index = cnt.index,
#                    columns = cnt.columns,
#                    )


coefs =  pd.read_csv(c_pth,sep = '\t',header = 0, index_col = 0)

inter = coefs.index.intersection(cnt.columns)

coefs = coefs.loc[inter]
cnt = cnt.loc[:,inter]

vals = np.dot(cnt.values,coefs.values).flatten()

crd = np.array([x.split('x') for x in cnt.index.values]).astype(float)


fig,ax = plt.subplots(1,1)
ax.scatter(crd[:,0],crd[:,1],s = 64, c = vals, cmap = plt.cm.Blues,marker = 'H')
ax.set_aspect("equal")
plt.show()




