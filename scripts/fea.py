#!/usr/bin/env python3

import numpy as np
import pandas as pd

import argparse as arp
import os.path as osp
import os
import sys

import matplotlib.pyplot as plt
from  matplotlib.colors import LinearSegmentedColormap as CMAP

from gprofiler import GProfiler

from typing import Tuple,Union,List

plt.rcParams['svg.fonttype'] = 'none'


""" FEA - Functional Enrichment Analysis

contains all functions to conduct and visualize
functional enrichment analysis using g:Profiler.
Also has a CLI-interface.

"""

def do_fea(genes : Union[List[str],np.ndarray],
           organism : str = "hsapiens",
           db : str = "GO:BP",
           )->pd.DataFrame:
    """ conduct FEA using g:Profiler"""

    gp = GProfiler(return_dataframe = True)
    res = gp.profile(organism = organism,
                     query = genes,
                     )

    if db.lower() != "all":
        if db in res['source'].values:
            gobps = res.loc[res['source'].values == db,:]
        else:
            print("[ERROR] : Not a valid DataBase")
            return pd.DataFrame([])

    gobps = gobps.loc[gobps['significant'].values,:]

    return gobps

def plot_fea_res(gobps : pd.DataFrame,
                 n_top : int = 20,
                 cmap : CMAP = plt.cm.YlOrRd,
                 )-> Tuple[plt.Figure,plt.Axes]:

    """visualize results from FEA"""

    def stretch(x,
                min_size = 0,
                max_size = 80,
                ):

        _range = x.max(keepdims = True) - x.min(keepdims =True)
        x_hat =  x - x.min(keepdims = True)
        x_hat /= _range

        _nrange = max_size - min_size

        out = x_hat*_nrange + min_size 

        return out



    label_dict = {"fontsize":20,
                  "family":"calibri"}
    
    standard_dict = {'family':'calibri'}
    n_top = np.min((n_top,gobps.shape[0]))

    y_pos = np.arange(n_top)[::-1] + 0.5
    p_vals = gobps['p_value'].values[0:n_top]
    x_pos = -np.log10(p_vals)

    fig,ax = plt.subplots(1,
                          1,
                          figsize = (13,8 /20 * n_top))

    cvals = stretch(x_pos,0,cmap.N) / 256
    cvals = cmap(cvals)

    delta = -np.diff(x_pos)[0] 
    ax.barh(y_pos,
            x_pos,
            color = cvals,edgecolor = 'black')

    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks(y_pos)
    ax.set_yticklabels(gobps['name'].values[0:n_top],**standard_dict)
    ax.set_aspect("equal")

    for k,pv in enumerate(x_pos):
        ax.text(pv + delta,y_pos[k],
                s = "{:0.3f}".format(pv),
                horizontalalignment = 'left',
                verticalalignment = 'center',
                **standard_dict,
                )

    for sp in ax.spines.values():
        sp.set_visible(False)

    ax.set_ylabel('Processes',**label_dict)
    ax.set_xlabel(r"$-\log10(p_{adj})$",**label_dict)

    return (fig,ax)


def main()->None:
    """Coefficient FEA CLI-interface"""

    prs = arp.ArgumentParser()
    aa = prs.add_argument

    aa("-cf","--coefficients",
       required = True,
       help = "path to coefficients",
       )
    aa("-od","--out_dir",
       default = None,
       help = "directory to save results",
       )

    aa("-nt","--n_top",
       default = 20,
       type = int,
       required = False,
       help = "number of top genes to print",
       )

    aa("-or","--organism",
       default = "hsapiens",
       help = "organism",
       )

    aa("-it","--image_type",
       choices = ["png",
                  "svg",
                  "jpg",
                  "gif",
                  ],
       default = "png",
       help = "organism",
       )


    aa("-cm","--color_map",
       default = "YlOrRd",
       type = str,
       help = "colormap; choose from" \
           " matplotlib's library",
       )


    aa("-db","--database",
       default = "GO:BP",
       required =False,
       type =str,
       help = "database to query against"\
       ". If non specified then GP:BP." \
       ' Specify "all" to query against all.',
       )

    args = prs.parse_args()

    if args.out_dir is None:
        out_dir = os.getcwd()
    else:
        out_dir = args.out_dir

    if not osp.exists(args.coefficients):
        print("[ERROR] : coefficient file does not exist. Exiting")
        sys.exit(-1)

    coef = pd.read_csv(args.coefficients,
                      sep = '\t',
                      header = 0,
                      index_col = 0)
    keep = pd.Index([x  for x in coef.index if x != "intercept"])

    coef = coef.loc[keep]

    genes = coef.index.values.tolist()

    gobps = do_fea(genes,
                   organism = args.organism,
                   db = args.database,
                   )
    bname = osp.basename(args.coefficients)
    bname = bname[0:min((20,len(bname)))]
    oname = osp.join(out_dir,"fea-" + bname)
    gobps.to_csv(oname + ".tsv",
                 sep = '\t',
                 header = True,
                 index =True,
                 )

    try:
        cmap = eval("plt.cm." + args.color_map)
        print(cmap)
    except:
        cmap = plt.cm.YlOrRd

    fig,ax = plot_fea_res(gobps,
                          n_top = args.n_top,
                          cmap = cmap,
                          )

    fig.savefig(oname + ".{}".format(args.image_type),
                transparent = True,
                )

if __name__ == "__main__":
    main()


