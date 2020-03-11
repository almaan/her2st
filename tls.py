#!/usr/bin/env python3

import pandas as pd
import numpy as np

import argparse as arp
import os.path as osp

import statsmodels.api as sm
from scipy.ndimage import gaussian_filter



def clean_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

def read_file(f,sep = '\t'):
    return pd.read_csv(f,header = 0,index_col = 0, sep = sep)

def joiner(pths,verbose = True):
    out = pd.DataFrame([])
    for k,pth in enumerate(pths):
        if verbose:
            print("[INFO] : Reading file {}".format(pth))
        out = pd.concat((out,read_file(pth)),sort = False)
        new_index = [str(k) + '_' + x for x in out.index.values]
        out.index = new_index
    out[out.isna()] = 0.0

    return out

def get_inflection_point(times : np.ndarray,
                         sigma : float = 10,
                         min_genes : int = 1,
                         n_fail : int = 100,
                         )-> int:

    f_times = gaussian_filter(times,sigma)
    f_d2 = gaussian_filter(np.gradient(np.gradient(f_times)),sigma)
    first = np.argmax(f_d2 > 0)
    f_d2[0:first] = 1
    ipoint = np.argmax(f_d2 <= 0)

    if ipoint < min_genes:
        ipoint = np.min((f_times,n_fail))
        print("[WARNING] : Did not find inflection point"
              " using {} genes instead.".format(ipoint))
    else:
        print("[INFO] : Using {} top genes".format(ipoint))

    return ipoint



prs = arp.ArgumentParser()
aa = prs.add_argument

aa("-c",
   "--count_files",
   nargs = '+')

aa("-p",
   "--proportion_files",
   nargs = '+')

aa("-o",
   "--out_dir",
   default = '/tmp')

aa("-n",
   "--n_top",
   default = 10,
   type = int)

args = prs.parse_args()

prop_pths = args.proportion_files
cnt_pths = args.count_files


#prop_pths = ["/home/alma/Documents/PhD/papers/HER2/res/major/ut_H{}_stdata_filtered.tsv/W.2020-03-04130406.862404.tsv".format(x) for x in range(1,4)]
#cnt_pths = ["/home/alma/Documents/PhD/papers/HER2/data/linnea-data/ut_H{}_stdata_filtered.tsv.gz".format(x) for x in range(1,4)]

prop_pths.sort()
cnt_pths.sort()

prop = joiner(prop_pths)
cnt = joiner(cnt_pths)


libSize = cnt.values.sum(axis = 1,keepdims= True)
nObs = cnt.values.sum(axis = 0, keepdims = True)

keep_genes = nObs.flatten() > 10
keep_spots = libSize.flatten() > 0

cnt = pd.DataFrame(np.divide(cnt.values, libSize,where = libSize > 0),
                   index = cnt.index,
                   columns = cnt.columns,
                   )

cnt = cnt.iloc[keep_spots,keep_genes]

inter = cnt.index.intersection(prop.index)

prop = prop.loc[inter,:]
cnt = cnt.loc[inter,:]




n_spots = cnt.shape[0]

type_1 = 'B-cells'
type_2 = 'T-cells'

pos_1 = np.argmax(prop.columns == type_1)
pos_2 = np.argmax(prop.columns == type_2)

# crd = np.array([x.split('x') for x \
#                 in  cnt.index.values]).astype(float)

jprod = np.zeros(n_spots)

for s in range(n_spots):

    vec = prop.values[s,:].reshape(-1,1)
    prod = np.dot(vec,vec.T)
    nprod = prod / prod.sum()

    jprod[s] = nprod[pos_1,pos_2]

jprod = pd.DataFrame(jprod,
                     index = cnt.index,
                     columns = ['probability'])

mod = sm.OLS(jprod,cnt)
res = mod.fit()

coefs = res.params
ordr = np.argsort(coefs.values)[::-1]
coefs = coefs.iloc[ordr]
pos = get_inflection_point(coefs.values)
coefs = coefs.iloc[0:pos]

coefs.to_csv(osp.join(args.out_dir,"tls-associated-2.tsv"),
             sep = '\t',
             header = True,
             index = True)

# marker_size = 5
# fig, ax = plt.subplots(1,11)

# ax[0].scatter(crd[:,0],
#             crd[:,1],
#             c = jprod.values.flatten(),
#             s = marker_size,
#               cmap = plt.cm.Reds)

# ax[0].set_aspect("equal")
# ax[0].set_title("Enrichement")

# for ii in range(1,11):

#     gene = top.index[ii-1]
#     vals = cnt[gene].values
#     vals = vals / vals.max()
#     rgba = np.zeros((n_spots,4))
#     rgba[:,2] = 1
#     rgba[:,3] = vals

#     ax[ii].scatter(crd[:,0],
#                 crd[:,1],
#                 c = jprod.values.flatten(),
#                 s = marker_size,
#                 alpha = 0.8,
#                 cmap = plt.cm.Reds)


#     ax[ii].scatter(crd[:,0],
#                 crd[:,1],
#                 c = rgba,
#                 s = marker_size,
#                 )



#     ax[ii].set_aspect("equal")
#     ax[ii].set_title(gene)

# for aa in ax:
#     clean_axes(aa)
# plt.show()

