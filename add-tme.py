#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os.path as osp
import os

import argparse as arp

from scipy.spatial.distance import cdist



def read_file(f,sep = '\t'):
    return pd.read_csv(f,header = 0,index_col = 0, sep = sep)

def get_crd(df):
    crd = np.zeros((df.shape[0],2))
    crd[:,0] = df['x'].values
    crd[:,1] = df['y'].values

    return crd


prs = arp.ArgumentParser()

prs.add_argument("-s",'--selections',nargs='+')

prs.add_argument("-o",'--out_dir',required = True)

prs.add_argument("-thr","--threshold",default = 2,type = float)

prs.add_argument("-x","--source",
                 default = ["invasive cancer","cancer in situ"],
                 nargs = '+')

prs.add_argument("-c","--col_name",default = 'label')


args = prs.parse_args()
# sel_dir =  "/tmp/spot_selections"
# sel_pths = list(filter(lambda f : f.split('.')[-1]  == "tsv", os.listdir(sel_dir)))
sel_pths = args.selections

patients = list(map(lambda x : osp.basename(x).split('_')[0],sel_pths))

out_dir = args.out_dir 
ind = args.source
col_name = args.col_name

mx = 0
tme_cutoff = args.threshold
eps = 0.1

ecs = []
fcs = []
crds = []

for p,sel_pth in enumerate( sel_pths ):

    sel = read_file(sel_pth)
    keep = np.any(pd.isna(sel),axis = 1) == False
    sel = sel.iloc[keep.values,:]
    crd = get_crd(sel)

    crds.append(crd)
    # is_tmr_idx = np.zeros(sel.shape[0])
    # is_non_idx = np.zeros(sel.shape[0])

    is_tmr_idx = np.array( [x in ind for x in sel[col_name].values ] )
    is_non_idx = np.array( [x not in ind for x in sel[col_name].values] )

    print(is_tmr_idx)

    is_tmr_name = sel.index.values[is_tmr_idx]
    is_non_name = sel.index.values[is_non_idx]

    crd_tmr  = crd[is_tmr_idx.flatten(),:]
    crd_non = crd[is_non_idx.flatten(),:]

    dmat = cdist(crd_non,crd_tmr)
    min_dist = np.min(dmat,axis = 1)

    if np.any(np.isnan(dmat.flatten())):
        print("[WARNING] : NaNs detected ")

    mx = np.max((mx,np.nanmax(min_dist)))

    sel['tumor_dist'] = np.zeros((crd.shape[0]))
    sel.loc[is_non_name,'tumor_dist'] = min_dist


    edgecolor = np.ones((crd.shape[0],4))
    edgecolor[:,3] = 1
    edgecolor[is_tmr_idx,0:3] = 0
    edgecolor[is_tmr_idx,0] = 1
    edgecolor[is_non_idx,0:3] = 0

    is_tme = (sel['tumor_dist'].values < tme_cutoff + eps) * is_non_idx

    edgecolor[is_tme,1] = 1

    ecs.append(edgecolor)
    fcs.append(sel['tumor_dist'].values)

    labels = np.array(['distal' for x in range( sel.shape[0] )],dtype = 'object')
    labels[is_tmr_idx] = 'cancer'
    labels[is_tme] = 'proximal'

    sel2 = sel.loc[:,:]

    sel2['label'] = labels

    sel2.to_csv(osp.join(out_dir, patients[p] + "-proximal-distal.tsv"),
                sep = '\t',
                header = True,
                index = True)


figsize = (5 * 3 + 0.5,5 * 3 + 0.5)
fig,ax = plt.subplots(3,3,figsize = figsize)
ax = ax.flatten()

for p in range(len(patients)):
    sc = ax[p].scatter(crds[p][:,0],
                       crds[p][:,1],
                       s = 60,
                       c = fcs[p].reshape(-1,),
                       edgecolor = ecs[p],
                       cmap = plt.cm.Blues,
                       vmin = 0,
                       vmax = mx,
                       )

    ax[p].set_aspect("equal")
    ax[p].set_title("Section : {}".format(patients[p]))
    ax[p].set_xticklabels([])
    ax[p].set_yticklabels([])


fig.colorbar(sc)
fig.tight_layout()

fig.savefig(osp.join(out_dir,
                     "tumor-dist-thrs-" + str(tme_cutoff) + ".png"))
plt.close("all")

