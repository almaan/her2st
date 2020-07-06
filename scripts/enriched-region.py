#!/usr/bin/env python3

import re
import os.path as osp
import argparse as arp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_1samp as ttest
from scipy.stats import wilcoxon
from scipy.stats import t as tdist

from typing import Union,Tuple

plt.rcParams['svg.fonttype'] = 'none'

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
    std = diff_dist.std(axis = 2)
    diff_dist = diff_dist.mean(axis =2)
    diff_dist /= std


    return {'vals':diff_dist,'regions':uni_labs}

def stretch(x,
            min_size = 0,
            max_size = 80,
            ):

    _range = x.max(keepdims = True) - x.min(keepdims =True)
    x_hat =  x - x.min(keepdims = True)
    x_hat /= _range

    print(x_hat.max(),x_hat.min())


    _nrange = max_size - min_size

    out = x_hat*_nrange + min_size 

    return out


def make_pallete(enr_res : dict,
                 feat_names : Union[np.ndarray,list],
                 min_size = 1,
                 max_size = 520,
                 )->Tuple[plt.Figure,plt.Axes]:


    xx = enr_res['vals']

    n_types = xx.shape[1]
    n_regions = xx.shape[0]
    figsize = (4 + n_regions * 0.5,
               2.5 + n_types * 0.5)

    fig, ax = plt.subplots(1,1,
                           figsize = figsize)

    x_pos = np.repeat(np.arange(n_regions),n_types) + 0.5
    y_pos = np.repeat(np.arange(n_types).reshape(1,-1),
                    n_regions,axis = 0).flatten()

    s = stretch(np.abs(xx.flatten()),min_size,max_size)
    is_up = (xx > 0).flatten()

    clr_up = np.array( [67, 191, 85] )
    clr_dw = np.array( [209, 56, 84] )

    rgba = np.zeros((s.shape[0],3))
    rgba[is_up,:] = clr_up / 256
    rgba[is_up == False,:]  = clr_dw / 256

    ax.scatter(x_pos,
            y_pos,
            s = s,
            c = rgba,
            )

    ax.set_xlabel("Regions")
    ax.set_ylabel("Types")
    ax.set_xticks(np.arange(n_regions)+0.5)
    ax.set_xlim([0,n_regions])
    ax.set_yticks(np.arange(n_types))
    ax.set_xticklabels(enr_res['regions'],
                       rotation = 90)
    ax.set_yticklabels(feat_names)

    ax.set_aspect("equal")

    fig.tight_layout()

    return(fig,ax)

def make_plots(diffs : np.ndarray,
               regions : np.ndarray,
               feats : np.ndarray,
               cmap = plt.cm.bone,
               )->dict:

    n_regions = diffs.shape[0]
    n_feats = diffs.shape[1]

    figsize = (n_regions*0.2,
               n_feats*1)

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
        image_type : str = "png",
        )->None:

    if tag is not None:
        tag = tag+"-"
    else:
        tag = ""

    if isinstance(labels,pd.DataFrame):
        labels = labels.values

    enr_res = enrichment(features.values,labels)

    fig, ax = make_pallete(enr_res,features.columns.values)


    fig.savefig(osp.join(out_dir,tag + "enrichment." + image_type),
                transparent = True)

def main(feats_pth : str,
         label_pth : str,
         label_col : str,
         out_dir : str,
         tag : str,
         subset_types : list,
         image_type : str,
         ):

    read_file = lambda f : pd.read_csv(f,
                                       index_col = 0,
                                       header = 0,
                                       sep = '\t')

    features = read_file(feats_pth)
    labels = read_file(label_pth)

    if subset_types is not None:
        keep_feat = [ x for x in features.columns.values if \
                    any([bool(re.search(z.lower(),x.lower())) for z in subset_types])]

        keep_feat = list(set(keep_feat))
        keep_feat.sort()
        keep_feat = pd.Index(keep_feat)
        features = features.loc[:,keep_feat]


    keep = labels.isna().any(axis=1).values == False
    labels = labels.iloc[keep,:]

    strint = lambda x : str(int(round(x)))

    new_idx = [ strint(x) + 'x' + strint(y) for \
                x,y in zip(labels['x'].values,labels['y'].values)]
    labels.index = pd.Index(new_idx)

    inter = features.index.intersection(labels.index)

    labels = labels.loc[inter,:]
    features = features.loc[inter,:]

    if label_col is not None:
        labels = labels[label_col].values

    run(features,labels,out_dir,tag,image_type)

if __name__ == "__main__":

    prs = arp.ArgumentParser()

    prs.add_argument("-f",
                     "--features")
    prs.add_argument("-l",
                     "--labels")
    prs.add_argument("-lc",
                     "--label_col",
                     default = None)
    prs.add_argument("-o",
                     "--output")
    prs.add_argument("-t",
                     "--tag",
                     default = None)
    prs.add_argument("-ext",
                     "--extension",
                     default = "png")
    prs.add_argument("-st",
                     "--subset_types",
                     default = None,
                     nargs = '+')


    args = prs.parse_args()

    main(args.features,
         args.labels,
         args.label_col,
         args.output,
         args.tag,
         args.subset_types,
         args.extension,
         )
