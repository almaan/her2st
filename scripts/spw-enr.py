#!/usr/bin/env python3

"""Spotwise enrichment analysis

Given a list of genes, e.g., representing
a set of genes associated to a pathway,
compute the spot-wise enrichment.

"""

import numpy as np
import pandas as pd

import os.path as osp
import os
import argparse as arp

from scipy.stats import fisher_exact

import matplotlib.pyplot as plt
from  matplotlib.colors import LinearSegmentedColormap as CMAP
from typing import Tuple


def read_file(f : str,
              sep : str = '\t',
              )->pd.DataFrame:

    return pd.read_csv(f,
                       header = 0,
                       index_col = 0,
                       sep = sep)


def get_crd(idx : pd.Index,
            )->np.ndarray:
    return np.array([x.split('x') for x in idx]).astype(np.float)


def get_enr(expr : pd.DataFrame,
            comp : pd.DataFrame,
            n_top : int = 100,
            )->np.ndarray:


    n_top = np.min((n_top,expr.shape[1]))
    n_spots = expr.shape[0]
    
    cnt_genes = expr.columns.values
    mu = expr.values.mean(axis = 0,keepdims = True)
    normalized_expr = expr.values - mu
    top_idx = np.argsort(normalized_expr,axis = 1)[:,::-1][:,0:n_top]
    all_set = set(cnt_genes)

    comp_set = set(comp)
    n_comp = len(comp_set)
    score = []

    for s in range(n_spots):
        cm = np.zeros((2,2))
        top = cnt_genes[top_idx[s,:]]
        top_set = set(top)
        inter = top_set.intersection(comp_set)
        union = top_set.union(comp_set)
        cm[0,0] = len(inter)
        cm[0,1] = len(comp_set.difference(top_set))
        cm[1,0] = len(top_set.difference(comp_set))
        cm[1,1] = len(all_set.difference(union))

        _,pvals = fisher_exact(cm)
        score.append( -np.log2(pvals))

    return np.array(score)

def clean_ax(axx : plt.Axes)->None:
    for sp in axx.spines.values():
        sp.set_visible(False)

    axx.set_xticks([])
    axx.set_xticklabels([])
    axx.set_yticks([])
    axx.set_yticklabels([])

def plot_spw_enr(crd : np.ndarray,
                 score : np.ndarray,
                 go_name : str,
                 marker_size : int = 300,
                 side_size : float = 10,
                 cmap : CMAP = plt.cm.Reds,
                 flip_y : bool = False,
                 )->Tuple[plt.Figure,plt.Axes]:
    figsize = (side_size,side_size)
    fig,ax = plt.subplots(1,
                          1,
                          figsize = figsize)

    vals = (score - score.min()) / (score.max() - score.min())

    ax.scatter(crd[:,0],
               crd[:,1],
               c = vals,
               s = marker_size,
               cmap = plt.cm.Reds,
               edgecolor = "gray",
               )

    ax.set_title(go_name)
    ax.set_aspect("equal")
    clean_ax(ax)

    if flip_y:
        ax.invert_yaxis()

    return fig,ax

def main()->None:
    prs = arp.ArgumentParser()
    aa = prs.add_argument

    aa("-c","--count_files",
       nargs = '+',
       required = True,
       help = "path to coefficients",
       )
    aa("-gs","--gene_set",
       required = True,
       help = "gene set to query against",
       )
    aa("-od","--out_dir",
       default = None,
       help = "directory to save results",
       )
    aa("-nt","--n_top",
       default = 100,
       type = int,
       required = False,
       help = "number of top genes use in analysis",
       )

    aa("-fy","--flip_y",
       default = False,
       required = False,
       action = "store_true",
       help = "number of top genes use in analysis",
       )


    aa("-ms","--marker_size",
       default = 300,
       type = int,
       required = False,
       help = "marker size",
       )

    aa("-ss","--side_size",
       default = 10,
       type = int,
       required = False,
       help = "figure side size",
       )


    aa("-cm","--color_map",
       default = "Reds",
       help = "colormap; choose from" \
           " matplotlib's library",
       )

    args = prs.parse_args()

    if args.out_dir is None:
        out_dir = os.getcwd()
    else:
        out_dir = args.out_dir

    sel_set_df = pd.read_csv(args.gene_set,
                          header = 0,
                          index_col = 0)


    go_name = osp.basename(args.gene_set)
    go_name = '.'.join(go_name.split('.')[0:-1])

    sel_col = ["name","gene"]
    for _col in sel_col:
        if _col in sel_set_df.columns:
            sel_set = sel_set_df[_col].values
            break 
        else:
            sel_set = None


    if sel_set is None:
        sel_set = sel_set_df.iloc[:,0].values

    for section in args.count_files:

        cnt = read_file(section)
        crd = get_crd(cnt.index)

        score = get_enr(cnt,sel_set,100)

        bname = osp.basename(section)
        bname = bname[0:min((20,len(bname)))]
        oname = osp.join(out_dir,"spwenr-" + go_name + "-" + bname)

        out_df = pd.DataFrame(dict(enrichment = score))
        out_df.index = cnt.index

        out_df.to_csv(oname + ".tsv",
                      sep = '\t',
                      header = True,
                      index = True,
                      )

        try:
            cmap = eval("plt.cm." + args.cmap)
        except:
            cmap = plt.cm.Reds

        fig,ax = plot_spw_enr(crd,
                              score,
                              go_name,
                              args.marker_size,
                              args.side_size,
                              cmap,
                              args.flip_y,
                              )

        fig.savefig(oname + ".png")

if __name__ == "__main__":
    main()
