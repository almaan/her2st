#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path as osp
import sys

from collections import OrderedDict
from operator import itemgetter

import argparse as arp

from typing import Tuple

def iprint(s : str) -> None:
    print("[INFO] : {}".format(s))

def epring(s : str) -> None:
    print("[ERROR] : {}".format(s))

def main():

    prs = arp.ArgumentParser()

    aa = prs.add_argument

    aa("-cf","--coef_file",
       required = True,
       type = str)
    aa("-nt","--number_top",
       required = False,
       default = None,
       help = "either a number or path" \
            " to the coef file only" \
            " containing the top ranked",
       )
    aa("-o","--out_dir",
       default ="/tmp")

    aa("-sg","--select_genes",
       nargs = "+",
       default = None,
       help = "list of genes to select" \
       " or path to file with genes",
       )

    aa("-t","--tag",
       required = False,
       default = None,
       )
    aa("-l","--lim",
       required = False,
       default = 0,
       type = float,
       )



    args = prs.parse_args()

    n_top = args.number_top

    coefs = pd.read_csv(args.coef_file,
                        sep = '\t',
                        index_col = 0,
                        header = 0)

    print(args.select_genes)

    if n_top is None:
        n_top = coefs.shape[0]
        iprint("will not mark top genes")
    elif osp.exists(n_top):
        iprint("reading file : {}".format(n_top))
        n_top = pd.read_csv(n_top,sep = '\t',
                            index_col = 0,
                            header = 0)
        inter = n_top.index.intersection(coefs.index)
        if inter.shape[0] == n_top.shape[0]:
            n_top = n_top.shape[0]
            iprint("number of top genes used : {}".format(n_top))
        else:
            eprint("coef filt and top file are not matched. Existing")
            sys.exit(-1)
    else:
        try:
            n_top = int(n_top)
        except:
            eprint("{} is not a path or number. Exiting.".format(n_top))
            sys.exis(-1)

    if args.select_genes is None:
        iprint("no genes selected for highlight")
    elif osp.exists(args.select_genes[0]):
        iprint("reading genes to highlight from : {}".format(args.select_genes))
        with open(args.select_genes[0],'r') as f:
            args.select_genes = [x.rstrip('\n').upper() for x in f.readlines()]
        iprint("genes prompted for highlight : {}".format(args.select_genes))
    elif isinstance(args.select_genes,list):
        args.select_genes = [x.upper() for x in args.select_genes]
        iprint("genes prompted for highlight : {}".format(args.select_genes))
    else:
        eprint("something went wrong in highlight selection. Exiting.")
        sys.exit(-1)

    fig,ax = make_rank_plot(coefs,
                            n_top,
                            args.select_genes,
                            args.lim,
                            )

    if args.tag is not None:
        tag = "-" + args.tag
    else:
        tag = ""

    out_pth = osp.join(args.out_dir,
                       "rank-plot" + tag + ".png"
                       )

    iprint("saving plot to : {}".format(out_pth))

    fig.savefig(out_pth)


def make_rank_plot(coefs : pd.DataFrame,
                   n_top : int,
                   sel_gen : list = None,
                   lim : float = 2,
                   )-> Tuple[plt.Figure,plt.Axes]:

    font_dict = {'fontfamily':'calibri',
                 'fontsize':20}


    srt = np.argsort(coefs.values.flatten())[::-1]
    coefs = coefs.iloc[srt,:]
    n_coefs = coefs.shape[0]

    color_top = ('red' if n_top < n_coefs else 'black')

    if sel_gen is not None:
        _coef_names = np.array(list(map(lambda x : x.upper(),coefs.index.values)))
        sel_rank = { x:np.argmax(x == _coef_names )  for \
                    x in sel_gen if x in _coef_names }

        sel_rank = OrderedDict(sorted(sel_rank.items(), key=itemgetter(1)))

        bg_rank = np.array([x for x in range(n_coefs) if\
                            x not in sel_rank.values()])
    else:
        sel_rank = {}
        bg_rank = np.arange(n_coefs)



    figsize = (15,15)

    fig,ax = plt.subplots(1,
                          1,
                          figsize = figsize)

    dy = 400*np.abs(np.diff(coefs.values.flatten()).mean())

    x_pos_bg = np.array(bg_rank)
    y_pos_bg = coefs.values.flatten()[bg_rank]

    ec = np.zeros((n_top,4))
    ec[:,0] = 1
    ec[:,3] = 0.01

    ax.scatter(x_pos_bg[0:n_top],
            y_pos_bg[0:n_top],
            color = color_top,
            s = 3,
            alpha = 0.4,
            )

    if n_top < n_coefs:
        ax.scatter(x_pos_bg[n_top::],
                y_pos_bg[n_top::],
                color = 'black',
                s = 1,
                marker = '.',
                alpha = 0.6,
                )
        
    if sel_rank is not None:
        y_pos_text = None
        crds = []
        for n,(k,v) in enumerate(sel_rank.items()):
            print(k,v)

            x_pos = v
            y_pos = coefs.values.flatten()[v]
            x_pos_text=  x_pos+1000 
            if y_pos_text is not None:
                x_pos_text=  x_pos+1400 
                y_high = y_pos + dy
                y_low =  y_pos - dy
                y_cand = np.linspace(y_low,y_high,100).tolist()
                get_dist = lambda x: np.abs(y_pos_text-x)
                dists = np.array(list(map(get_dist,y_cand)))
                if np.all(dists/dy > lim):
                    y_pos_text = y_high
                elif np.any(dists /dy > lim ):
                    y_pos_text = y_cand[np.argmax(np.where(dists/dy > lim)[0])]
                else:
                    y_pos_text = y_low
                    # y_pos_text = y_low

            else:
                y_pos_text = dy + y_pos 

            crds.append([x_pos_text,y_pos_text])

            ax.plot((x_pos_text,x_pos),
                    (y_pos_text,y_pos),
                    color = 'black'
                    )

            ax.scatter(x_pos,
                    y_pos,
                    s = 60,
                    c = 'red',
                    edgecolor = 'black',
                    alpha = 0.8,
                    )
            ax.text(s = "{} : {}".format(k,v+1),
                    x = x_pos_text,
                    y = y_pos_text,
                    **font_dict,
                    )

    ax.tick_params(axis='both',
                which='major',
                labelsize=12)
    ax.tick_params(axis ='x',
                which = 'major',
                rotation = 45)

    ax.set_xlabel("Rank",
                  **font_dict)
    ax.set_ylabel("Coefficient Value",
                **font_dict)

    mn = coefs.values.flatten().min()
    mx = coefs.values.flatten().max()
    ax.set_ylim([1.2**(-np.sign(mn))*mn,1.2**np.sign(mx)*mx])


    for sp in ['top','right']:
        ax.spines[sp].set_visible(False)

    for sp in ['left','bottom']:
        ax.spines[sp].set_lw(1.1)


    return fig,ax

if __name__ == "__main__":
    main()
