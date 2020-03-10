#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

import os.path as osp
import re
import argparse as arp


def read_file(f,sep = '\t'):
    return pd.read_csv(f,header = 0,index_col = 0, sep = sep)

def joiner(pths):
    out = pd.DataFrame([])
    for k,pth in enumerate(pths):
        out = pd.concat((out,read_file(pth)),sort = False)
        new_index = [str(k) + '_' + x for x in out.index.values]
        out.index = new_index
    out[out.isna()] = 0.0

    return out

def clean_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])


prs = arp.ArgumentParser()
aa = prs.add_argument

aa("-c",
   "--count_files",
   nargs = '+')

aa("-p",
   "--proportion_files",
   nargs = '+')

aa("-t",
   "--types",
   nargs = "+",
   default = None)

aa("-o",
   "--out_dir",
   default = '/tmp')

aa("-n",
   "--n_top",
   default = 10,
   type = int)

args = prs.parse_args()
# prop_pths = ["/home/alma/Documents/PhD/papers/HER2/res/subset/ut_H2_stdata_filtered.tsv/W.2020-03-06004614.163337.tsv"]
# cnt_pths = ["/home/alma/Documents/PhD/papers/HER2/data/linnea-data/ut_H2_stdata_filtered.tsv.gz"]

prop_pths = args.proportion_files
cnt_pths = args.count_files

n_top = args.n_top

prop_pths.sort()
cnt_pths.sort()

cnt = joiner(cnt_pths)
prop = joiner(prop_pths)

if args.types is not None:
    sel_cols = [ x for x in prop.columns.values if \
                 any([bool(re.search(z.lower(),x.lower())) for z in args.types])]

    sel_cols = pd.Index(sel_cols)
    prop = prop.loc[:,sel_cols]

inter = cnt.index.intersection(prop.index)
cnt = cnt.loc[inter,:]

keep1 = (cnt.values.sum(axis = 0) > 20)

cnt = cnt.iloc[:,keep1]

prop = prop.loc[inter,:]

crd = np.array([x.split('_')[1].split('x') for x \
                in  cnt.index.values]).astype(float)

cell_types = prop.columns.values
top_genes = dict()

for tp in cell_types:
    corr = cdist(cnt.values.T,
                 prop[tp].values.reshape(1,-1),
                 metric = "correlation")

    corr =  - corr + 1

    corr = corr.flatten()
    corr[np.isnan(corr)] = np.nanmin(corr)

    srt = np.argsort(corr)
    srt = srt.flatten()[::-1][0:1000]
    top = cnt.columns.values[srt]

    top_genes.update({tp:top})


marker_size = 40
n_spots = cnt.shape[0]
n_top = args.n_top
n_types = cell_types.shape[0]
side_length = 3

plt_args = dict(s = marker_size,
                edgecolor = 'none',
                marker = 'o',
                )

for tp in cell_types:
    fig,ax = plt.subplots(3,n_top,figsize = (n_types * side_length, side_length * 3))

    ax[0,0].set_ylabel("Gene expression")
    ax[1,0].set_ylabel("Type proportions")
    ax[2,0].set_ylabel("Joint Overlay")

    fig.suptitle(tp)

    for k in range(n_top):
        gene = top_genes[tp][k]
        cnt_clr = np.zeros((n_spots,4))
        prop_clr = np.zeros(cnt_clr.shape)

        cnt_clr[:,2] = 1
        cnt_clr[:,3] = cnt[gene].values
        cnt_clr[:,3] /= cnt_clr[:,3].max()
        prop_clr[:,0] = 1
        prop_clr[:,3] = prop[tp].values / prop[tp].values.max()

        joint_clr_1 = prop_clr[:,:]
        joint_clr_1[:,3] /= 2
        joint_clr_2 = cnt_clr[:,:]
        joint_clr_2[:,3] /= 2

        ax[0,k].scatter(crd[:,0],
                        crd[:,1],
                        c = cnt_clr,
                        **plt_args,
                        )

        ax[0,k].set_title(gene)

        ax[1,k].scatter(crd[:,0],
                        crd[:,1],
                        c = prop_clr,
                        **plt_args,
                        )

        # ax[1,k].set_title(tp.replace("_","\n"),
                          # fontsize = 15)

        ax[2,k].scatter(crd[:,0],
                        crd[:,1],
                        c = joint_clr_1,
                        **plt_args,
                        )

        ax[2,k].scatter(crd[:,0],
                        crd[:,1],
                        c = joint_clr_2,
                        **plt_args)

        ax[0,k].set_aspect("equal")
        ax[1,k].set_aspect("equal")
        ax[2,k].set_aspect("equal")

    for aa in ax.flatten(): clean_axes(aa)

    fig.savefig(osp.join(args.out_dir,
                         tp.replace('_','-').replace(' ','-') + '.png'))


