import numpy as np
import pandas as pd

import argparse as arp
import os.path as osp

import matplotlib.pyplot as plt

from gprofiler import GProfiler

from typing import Tuple

import matplotlib

CMAP = matplotlib.colors.LinearSegmentedColormap

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



def do_fea(genes : list,
           organism : str = "hsapiens",
           db : str = "GO:BP",
           )->pd.DataFrame:

    gp = GProfiler(return_dataframe = True)
    res = gp.profile(organism = organism,
                     query = genes,
                     )

    gobps = res.loc[res['source'].values == db,:]
    gobps = gobps.loc[gobps['significant'].values,:]

    return gopbs

def plot_fea_res(gobps : pd.DataFrame,
                 n_top : int = 20,
                 cmap : CMAP = plt.cm.YlOrRd,
                 )-> Tuple[plt.Figure,plt.Axes]:

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


# Bean plots are not available in
# matplotlib, this is a
# custom implementation utilzing
# the violinplot function

def bean_plot(ax : plt.Axes,
              vals : list,
              cmap : CMAP = plt.cm.Set1,
              edgecolor : str = "black",
              max_width = 0.15,
              )->plt.Axes:

    xs = np.arange(len(vals))
    vps = ax.violinplot(vals,
                        showextrema=False,
                        positions = xs,
                        )

    ax.set_xticks(xs)

    for k,patch in enumerate(vps["bodies"]):
        patch.set_facecolor(cmap(k))
        patch.set_edgecolor(edgecolor)
        patch.set_alpha(1.0)


    for center,dist in zip(xs,vals):
        points,counts = np.unique(dist,return_counts = True)
        delta = max_width / 2

        _xp = (center - delta,
               center + delta)

        _xp = np.array(_xp)

        mu = np.median(dist)

        ax.plot((center - delta * 2,
                 center + delta * 2),
                (mu,mu),
                linewidth = 2,
                color  = "red",
                )

        ax.scatter(center,
                   mu,
                   facecolor = "red",
                   edgecolor = "black",
                   zorder = np.inf,
        )

        for p in points:

            _yp = (p,p)

            ax.plot(_xp,
                    _yp,
                    color = "black",
                    linewidth=1.5,
                    alpha =0.3,
                    )

    return vps

