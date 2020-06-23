#!/usr/bin/env python3

import pandas as pd
import numpy as np

import yaml

import os.path as osp
import argparse as arp

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from typing import Union,Dict,Tuple,List

plt.rcParams['svg.fonttype'] = 'none'

def make_plots(features_path : str,
               labels_paths : Dict[str,str],
               out_dir : str,
               cmap : Dict[int,Tuple[int,int,int]],
               labels_columns : Dict[str,str],
               select_types : Dict[str,str],
               annotate : bool,
               )-> None:

    read_file = lambda f : pd.read_csv(f,
                                       sep = '\t',
                                       header = 0,
                                       index_col = 0)


    for lab_type,l_pth in labels_paths.items():

        feats = read_file(features_path)
        labs = read_file(l_pth)

        keep_labs = (pd.isna(labs)\
                     .any(axis =1)\
                     .values == False)

        labs = labs.iloc[keep_labs,:]

        inter = feats.index.intersection(labs.index)

        if len(inter) < 1:
            fmt = lambda x: str(int(round(x)))
            new_idx = pd.Index([fmt(x) + \
                                "x" + \
                                fmt(y) for \
                                x,y in zip(labs["x"].values,
                                           labs["y"].values)])

            labs.index = new_idx

        inter = feats\
               .index\
               .intersection(labs.index)

        feats = feats.loc[inter,:]
        labs = labs.loc[inter,:]

        crd = labs[["pixel_x",
                    "pixel_y"]].values

        lab_col = labels_columns[lab_type]

        for sel_type,display_name in select_types.items():
            edgecolor = np.zeros((feats.shape[0],4))
            patches = list()

            fig,ax = plt.subplots(1,1,figsize = (12,12))

            for lab in cmap[lab_type].keys():
                pos = labs[lab_col].values == lab
                edgecolor[pos,0:3] = cmap[lab_type][lab]
                edgecolor[pos,3] = 255.0 * 0.8

                _clr = np.array(cmap[lab_type][lab] ) / 255.0
                _patch = mpatches.Patch(color = _clr,
                                        label = f"{lab_type} : {lab}")
                patches.append(_patch)

            edgecolor /= 255.0

            scale = lambda x : x / x.max() * 0.8

            ax.scatter(x = crd[:,0],
                    y = crd[:,1],
                    c = None,
                    s = 360,
                    edgecolor = "gray",
                    linewidth = 2,
                    marker = "o",
                    )

            ax.scatter(x = crd[:,0],
                    y = crd[:,1],
                    c = scale(feats[sel_type].values),
                    s = 330,
                    cmap = plt.cm.gray_r,
                    edgecolor = edgecolor,
                    linewidth = 3,
                    marker = "o",
                    vmin = 0,
                    vmax = 1,
                    )

            ax.set_aspect("equal")
            if lab_type == "cluster":
                ax.invert_yaxis()

            if annotate : ax.legend(handles = patches)

            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_yticklabels([])
            ax.set_title(f"{display_name}")

            for sp in ax.spines.values():
                sp.set_visible(False)

            fig.savefig(osp.join(out_dir,"{}-{}-based-spotplot.png"\
                        .format(sel_type,lab_type)),
                        transparent = True,
                        )
            plt.close("all")

    return None


def main():

    prs = arp.ArgumentParser()
    aa = prs.add_argument

    aa("-c","--config_file",
       type = str,
       required = True,
       )

    aa("-o","--output_dir")

    args = prs.parse_args()

    with open(args.config_file,"r+") as f:
        config = yaml.load(f)

    make_plots(features_path = config["data"]["features_path"],
               labels_paths = config["data"]["labels_paths"],
               out_dir = args.output_dir,
               cmap = config["cmap"],
               select_types = config["select_types"],
               annotate = config["annotate"],
               labels_columns = config["labels_columns"]
               )


if __name__ == "__main__":
    main()
