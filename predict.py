#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os.path as osp
import argparse as arp

def iprint( s : str,
            ) -> None:
    print("[INFO] : {}".format(s))

def eprint( s : str,
            ) -> None:
    print("[ERROR] : {}".format(s))



prs = arp.ArgumentParser()

aa = prs.add_argument

aa("-cn","--count_files",
   required = True,
   nargs = '+',
   help = 'count files',
   )

aa('-cf',"--coefficients",
   required = True,
   help = 'coefficients path',
   type = str,
   )

aa("-eps","--threshold",
   default = None,
   type = float,
   help = 'threshold for prediction values',
   )

aa("-o","--out_dir",
   default = "/tmp/",
   help = 'where to save output',
   )

aa("-cm","--color_map",
   default = "PuRd",
   type = str,
   )

aa("-ms","--marker_size",
   default = 64,
   type = int,
   )

aa("-mt","--marker_type",
   default = "H",
   type = str,
   )

args = prs.parse_args()

try:
    args.color_map = eval("plt.cm." + args.color_map)
except:
    eprint("{} is not a color map. Using default".format(args.color_map))
    args.color_map = plt.cm.PuRd

for pth in args.count_files:

    coefs =  pd.read_csv(args.coefficients,
                         sep = '\t',
                         header = 0,
                         index_col = 0)

    cnt = pd.read_csv(pth,
                       sep = '\t',
                       header = 0,
                       index_col = 0)

    cnt = pd.DataFrame(cnt.values / cnt.values.sum(axis = 1,
                                                   keepdims = True),
                    index = cnt.index,
                    columns = cnt.columns,
                    )

    if 'intercept' in coefs.index:
        iprint("intercept included in prediction")
        cnt['intercept'] = np.ones(cnt.shape[0])

    inter = coefs.index.intersection(cnt.columns)

    iprint("{} / {} genes present in count data".format(inter.shape[0],
                                                        cnt.shape[1],
                                                        ))

    coefs = coefs.loc[inter]
    cnt = cnt.loc[:,inter]

    vals = np.dot(cnt.values,coefs.values)
    vals = vals.flatten()

    if args.threshold is not None:
        iprint("using cutoff : {}".format(args.threshold))
        vals[vals < args.threshold] = 0.0
    else:
        _threshold = "None"

    edgecolor = np.zeros((vals.shape[0],4))
    edgecolor[:,3] = 0.4

    crd = np.array([x.split('x') for \
                    x in cnt.index.values]).astype(float)

    figsize = (10,10)
    fig,ax = plt.subplots(1,1,
                          figsize = figsize)

    ax.scatter(crd[:,0],
            crd[:,1],
            s = args.marker_size,
            c = vals,
            marker = args.marker_type,
            cmap = args.color_map,
            edgecolor = edgecolor,
            )

    ax.set_aspect("equal")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    bname = osp.basename(pth).replace(".tsv","").replace(".gz","")

    if len(bname) >= 20:
        bname = bname[0:20]

    out_pth = osp.join(args.out_dir,
                         "tls-pred-cutoff-{}-{}.png".format(str(_threshold),bname))

    iprint("saving results of {} to file >> {}".format(pth,out_pth))
    fig.savefig(out_pth)
    plt.close("all")


