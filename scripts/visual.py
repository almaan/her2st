#!/usr/bin/env python3

import pandas as pd
import numpy as np

import os.path as osp
import os
import sys

import argparse as arp

import matplotlib.pyplot as plt

plt.close("all")

def read_file(f,sep = '\t')->pd.DataFrame:
    return pd.read_csv(f,header = 0,index_col = 0, sep = sep)


def _normalize(x : np.ndarray)->np.ndarray:
    # as of now no normalization
    return x
    # return np.log2(x+1)

def get_crd(idx : pd.Index):
    return np.array([x.split('x') for x in idx]).astype(np.float)

def clean_ax(axx : plt.Axes)->None:
    for sp in axx.spines.values():
        sp.set_visible(False)

    axx.set_xticks([])
    axx.set_xticklabels([])
    axx.set_yticks([])
    axx.set_yticklabels([])

def hex2rgb(H : str)->np.ndarray:
    try:
        h = H.lstrip("#").upper()
        return np.array(tuple(int(h[i:i+2], 16) for i in (0, 2, 4))) / 256
    except:
        print("[EROR] : not a proper hex-value. Defualt color will be used")
        return np.array([235, 186, 52]) / 256



prs = arp.ArgumentParser()

aa = prs.add_argument

aa("-c","--count_file",
   required = True,
   type = str,
   help = "full path to count file",
   )

aa("-o","--out_dir",
   required = False,
   type = str,
   default = os.getcwd(),
   help = ' '.join(("Directory to save files to.",
           "CWD will be used if none specified")),
   )

aa("-ms","--marker_size",
   default = 30,
   type = int,
   help = "marker size"
   )

aa("-g","--genes",
   nargs = "+",
   help = "list of genes to plot. Or path to list."
   )

aa("-mc","--marker_color",
   default = None,
   help = "must be hex value"
   )

aa("-ss","--side_size",
   default = 10,
   type = int,
   help = "size of image"
   )

aa("-fs","--font_size",
   default = 40,
   type = int,
   help = "font size"
   )

aa("-tb","--transparent_background",
   action = 'store_true',
   default = False,
   help = "include flag for transparent background",
   )



args = prs.parse_args()

if not osp.isdir(args.out_dir):
    print("[INFO] : creating directory >> {}".format(args.out_dir))
    os.mkdir(args.out_dir)

cnt = read_file(args.count_file)
crd = get_crd(cnt.index)

if args.marker_color is not None:
    color =  hex2rgb(args.marker_color)
else:
    color = np.array([235, 186, 52]) / 256

marker_size = args.marker_size

if osp.exists(args.genes[0]):
    try:
        with open(args.genes[0],"r") as f:
            sel_genes = f.readlines()
            sel_genes = [x.rstrip('\n').upper() for x in sel_genes]
    except:
        print("[ERROR] : could not read gene file.Exiting.")
        sys.exit(-1)
else:
    sel_genes = args.genes

cnt.columns = pd.Index([x.upper() for x in cnt.columns])

for gene in sel_genes: 

    if gene not in cnt.columns.values:
        print("[INFO] : {} not present. Will skip.".format(gene))
        continue
    else:
        print("[INFO] : Rendering {}".format(gene))

    figsize = (args.side_size,args.side_size)
    fig,ax = plt.subplots(1,1,figsize = figsize)

    rgb = np.zeros((cnt.shape[0],4))
    rgb[:,0:3] = color
    n_vals = _normalize(cnt[gene].values).reshape(-1,)
    rgb[:,3] = n_vals / n_vals.max()

    ax.scatter(crd[:,0],
               crd[:,1],
               c = rgb,
               s = args.marker_size * 10,
               edgecolor = "gray",
               )

    ax.set_title(gene,
                 fontsize = args.font_size)
    ax.set_aspect("equal")
    clean_ax(ax)

    bname = osp.basename(args.count_file)
    out_pth = osp.join(args.out_dir,
                       bname + "-viz-{}.png".format(gene))
    fig.savefig(out_pth,
                transparent = args.transparent_background)

