#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import argparse as arp
import os.path as osp

def read_file(f,sep = '\t'):
    return pd.read_csv(f,header = 0,index_col = 0, sep = sep)


def run(prop_pths : str,
        out_pth : str,
        marker_size : float = 80,
        edgecolor : str = "black",
        transparent_background  : bool = False,
        hide_title : bool = False,
        ):


    prop = read_file(prop_pths)
    n_spots = prop.shape[0]

    type_1 = 'B-cells'
    type_2 = 'T-cells'

    pos_1 = np.argmax(prop.columns == type_1)
    pos_2 = np.argmax(prop.columns == type_2)

    crd = np.array([x.split('x') for x \
                    in  prop.index.values]).astype(float)

    jprod = np.zeros(n_spots)


    jprod = np.zeros(n_spots)

    for s in range(n_spots):
        vec = prop.values[s,:].reshape(-1,1)
        prod = np.dot(vec,vec.T)
        nprod = prod / prod.sum()
        N = prod.shape[0]

        jprod[s] = nprod[pos_1,pos_2] * 2
        jprod[s] -= (nprod.sum() / (0.5*(N**2 + N)))


    value_list = [prop[type_1].values,prop[type_2].values,jprod.flatten()]

    max_val = 0
    for ii in value_list:
        mx = ii.max()
        if mx > max_val:
            max_val = mx

    color_list = [[163, 44, 46],
                [44, 163, 78],
                [120, 44, 122],
                ]

    color_list = np.array(color_list)
    titles = [type_1,type_2,"TLS-score"]

    figsize = (42,14)
    fig,ax = plt.subplots(1,3,figsize = figsize)

    for k,_vals in enumerate(value_list):

        vals = _vals.flatten()
        vals -= vals.min()
        rgba = np.zeros((vals.shape[0],4))
        rgba[:,0:3] = (color_list[k,:] / 256)
        rgba[:,3] = vals / vals.max()

        edgecolor = np.zeros(rgba.shape)
        edgecolor[:,3] = 0.4

        ax[k].scatter(crd[:,0],
                    crd[:,1],
                    c = rgba,
                    s = marker_size,
                    edgecolor = edgecolor,
                    )

        ax[k].set_aspect("equal")
        if not hide_title : ax[k].set_title(titles[k])

    for aa in ax:
        aa.invert_yaxis()
        aa.set_xticklabels([])
        aa.set_yticklabels([])
        aa.set_xticks([])
        aa.set_yticks([])
        for sp in aa.spines.values():
            sp.set_visible(False)

    fig.savefig(out_pth,
                transparent = transparent_background,
                )

def main():

    prs = arp.ArgumentParser()
    aa = prs.add_argument

    aa("-p",
    "--proportion_files",
       )

    aa("-o",
    "--out_dir",
    default = '/tmp')

    aa("-t",
    "--tag",
    default = None,
    )

    aa("-it",
       "--image_type",
       type = str,
       choices = ["png",
                  "svg",
                  "jpg",
                  ],
       default = "png",
       help = "",
       )

    aa("-tbg",
        "--transparent_background",
        default = False,
        action = "store_true",
       )

    aa("-ht",
        "--hide_title",
        default = False,
        action = "store_true",
       )


    aa("-ms",
       "--marker_size",
       type = float,
       default = 280,
       )


    args = prs.parse_args()

    if args.tag is not None:
        args.tag = args.tag + "-"
    else:
        args.tag = ""


    out_pth = osp.join(args.out_dir,
                       args.tag +\
                       "tls-visualize." + \
                       args.image_type)

    run(args.proportion_files,
        out_pth,
        marker_size = args.marker_size,
        transparent_background = args.transparent_background,
        hide_title = args.hide_title,
        )

if __name__ == "__main__":
    main()

