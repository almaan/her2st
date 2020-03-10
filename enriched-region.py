#!/usr/bin/env python3


import os.path as osp
import argparse as arp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_1samp as ttest
from scipy.stats import wilcoxon
from scipy.stats import t as tdist

from typing import Union

# def gene_1(x,y,mu = 10):
#     return mu + x + np.random.normal(0,1,x.shape)

# def gene_2(x,y,mu = 10):
#     return mu + y + np.random.normal(0,1,y.shape)

# def gene_3(x,y,mu = 10):
#     return mu + (x-15)**2 + (y-15)**2 + np.random.normal(0,1,x.shape)

# def quadrants(x,y):
#     labs = np.zeros(x.shape[0]).astype(int)
#     for ii in range(x.shape[0]):
#         if x[ii] < 15:
#             if y[ii] < 15:
#                 labs[ii] = 1
#             else:
#                 labs[ii] = 2
#         else:
#             if y[ii] < 15:
#                 labs[ii] = 3
#             else:
#                 labs[ii] = 4
#     return labs


def _generate_null(xmat,
                   labels,
                   n_shuffle,
                   )->np.ndarray:

    uni_labs = np.unique(labels)
    mu_null = np.zeros((uni_labs.shape[0],
                        xmat.shape[1],
                        n_shuffle,
                        ))

    for it in range(n_shuffle):
        idx = np.random.permutation(xmat.shape[0])
        for num,lab in enumerate(uni_labs):
            mu_null[num,:,it] = xmat[idx[labels==lab],:].mean(axis = 0)

        if it % 1000 == 0 and it >= 1000:
            print("\r Iter : {}".format(it),end="")

    print("\n",end="")

    return mu_null

def enrichment(xmat : np.ndarray,
               labels : np.ndarray,
               n_shuffle : int = 10000,
               )->dict:

    uni_labs = np.unique(labels,)
    n_regions = uni_labs.shape[0]
    n_feature = xmat.shape[1]
    true_mean = np.zeros((n_regions,n_feature,1))

    for num,lab in enumerate(uni_labs):
        true_mean[num,:,0] = xmat[labels == lab,:].mean(axis=0)

    null_dist = _generate_null(xmat,labels,n_shuffle)

    diff_dist = true_mean - null_dist

    return {"diffs":diff_dist,"regions":uni_labs}

def do_test(diffs : np.ndarray,
            alt_mu : int = 0,
            ax : int = 2,
            )-> np.ndarray:

    return ttest(diffs,
                 popmean = alt_mu,
                 axis = ax,
                 )

def get_stats(diffs : np.ndarray,
              alpha = 0.05,
              ):

    stds = diffs.std(axis =2)
    mus = diffs.mean(axis =2)
    n = diffs.shape[2]
    tstar = tdist.ppf((1-alpha/2),n -1)

    stats = dict(mus = mus,
                 errs = tstar*stds/n**0.5,
                 )

    return stats

def make_plots(diffs : np.ndarray,
               regions : np.ndarray,
               feats : np.ndarray,
               cmap = plt.cm.bone,
               )->dict:

    n_regions = diffs.shape[0]
    n_feats = diffs.shape[1]

    figsize = (n_regions*3,15)

    vizlist = dict()

    mx = np.abs(diffs).max()

    for k,region in enumerate(regions):

        fig,ax = plt.subplots(1,1,figsize = figsize)
        sub_diff_dist = diffs[k,:,:].T

        bp = ax.boxplot(sub_diff_dist,
                        sym = "",
                        patch_artist = True,
                        )


        for ft in range(n_feats):
            ypos = bp['caps'][ft*2].get_ydata()[0]
            ypos -= mx*0.01

            ax.text(x = ft +1,
                    y = ypos,
                    s = feats[ft],
                    fontdict = {"rotation":90,
                                "horizontalalignment":"center",
                                "verticalalignment":"top",
                                },
                    )

        clrs = np.linspace(0,cmap.N,n_feats) / cmap.N

        for c,patch in enumerate(bp['boxes']):
            patch.set_facecolor(cmap(clrs[c]))


        # ax.set_xticks(np.arange(n_feats) + 1)
        ax.set_xticks([])
        ax.set_xticklabels([])

        # ax.set_xticklabels(feats,
        #                    rotation = 90)
        ax.axhline(y = 0,
                   linestyle = 'dashed',
                   linewidth = 1,
                   color = 'black',
                   )

        ax.set_title("Region : {}".format(region))
        ax.set_ylabel(r"$\Delta$")
        ax.set_ylim([-mx,mx])

        for sp in ax.spines.values():
            sp.set_visible(False)

        fig.tight_layout()

        vizlist.update({region:(fig,ax)})

    return vizlist

def run(features : pd.DataFrame,
        labels : Union[np.ndarray,pd.DataFrame],
        out_dir : str = "/tmp",
        tag : str = None,
        )->None:

    if tag is not None:
        tag = tag+"-"
    else:
        tag = ""

    if isinstance(labels,pd.DataFrame):
        labels = labels.values

    enr_res = enrichment(features.values,labels)

    vizlist = make_plots(enr_res['diffs'],
                        enr_res['regions'],
                        feats = features.columns.values)

    for key,val in vizlist.items():
        val[0].savefig(osp.join(out_dir,tag + "region-" + key + ".png"))

def main(feats_pth : str,
         label_pth : str,
         label_col : str,
         out_dir : str,
         tag : str,
         ):

    read_file = lambda f : pd.read_csv(f,
                                       index_col = 0,
                                       header = 0,
                                       sep = '\t')

    features = read_file(feats_pth)
    labels = read_file(label_pth)
    keep = labels.isna().any(axis=1).values == False
    labels = labels.iloc[keep,:]

    strint = lambda x : str(int(round(x)))

    new_idx = [ strint(x) + 'x' + strint(y) for x,y in zip(labels['x'].values,labels['y'].values)]
    labels.index = pd.Index(new_idx)

    inter = features.index.intersection(labels.index)

    labels = labels.loc[inter,:]
    features = features.loc[inter,:]

    if label_col is not None:
        labels = labels[label_col].values

    run(features,labels,out_dir,tag)

if __name__ == "__main__":

    prs = arp.ArgumentParser()

    prs.add_argument("-f","--features")
    prs.add_argument("-l","--labels")
    prs.add_argument("-lc","--label_col",default = None)
    prs.add_argument("-o","--output")
    prs.add_argument("-t","--tag",default = None)

    args = prs.parse_args()

    main(args.features,args.labels,args.label_col,args.output,args.tag)
