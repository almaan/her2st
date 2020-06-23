#!/usr/bin/env python3

"""TLS-score prediction"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os.path as osp
import argparse as arp
import sys

def iprint( s : str,
            ) -> None:
    print("[INFO] : {}".format(s))

def eprint( s : str,
            ) -> None:
    print("[ERROR] : {}".format(s))


def main()->None:
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

    aa("-tp","--transpose",
    default = False,
    action = 'store_true',
    )

    aa("-sn","--split_names",
    nargs = '+',
    default = None,
    )

    args = prs.parse_args()

    try:
        args.color_map = eval("plt.cm." + args.color_map)
    except:
        eprint("{} is not a color map. Using default".format(args.color_map))
        args.color_map = plt.cm.PuRd


    maxV = -np.inf
    minV = np.inf
    valS = []
    crdS = []

    for pth in args.count_files:
        try:
            coefs =  pd.read_csv(args.coefficients,
                                sep = '\t',
                                header = 0,
                                index_col = 0)

            cnt = pd.read_csv(pth,
                            sep = '\t',
                            header = 0,
                            index_col = 0)

            if args.transpose:
                cnt = cnt.T

            if args.split_names is not None:
                new_cols = pd.Index([x.split(args.split_names[0])[int(args.split_names[1])]for\
                                    x in cnt.columns])
                cnt.columns = new_cols



            cnt = pd.DataFrame(cnt.values / cnt.values.sum(axis = 1,
                                                        keepdims = True),
                            index = cnt.index,
                            columns = cnt.columns,
                            )

            std = cnt.values.std(axis=0,keepdims = True)

            cnt = pd.DataFrame(np.divide(cnt.values,std ,where = std > 0),
                    index = cnt.index,
                    columns = cnt.columns,
                    )


            if 'intercept' in coefs.index:
                iprint("intercept included in prediction")
                cnt['intercept'] = np.ones(cnt.shape[0])

            cnt.columns = pd.Index([x.upper().replace(" ",'') for x in cnt.columns.values])
            coefs.index = pd.Index([x.upper().replace(" ","") for x in coefs.index.values])

            uniq_genes,gene_occ = np.unique(cnt.columns.values,return_counts =True)
            cnt = cnt.loc[:,uniq_genes[gene_occ == 1]]

            inter = pd.Index(set(coefs.index.intersection(cnt.columns)))

            iprint("{} / {} genes present in count data".format(inter.shape[0],
                                                                coefs.shape[0],
                                                                ))
            coefs = coefs.loc[inter,:]


            cnt = cnt.loc[:,inter]
            uu,cc = np.unique(cnt.columns.values,return_counts = True)
            diff = coefs.index.difference(cnt.columns)

            val = np.dot(cnt.values,coefs.values)
            val = val.flatten()
            iprint("Dynamic range of values : MIN({:0.4f}) , MAX({:0.4f})".format(val.min(),val.max()))

            crd = np.array([x.replace('X',"").split('x') for \
                            x in cnt.index.values]).astype(float)

            maxV = np.max((maxV,val.max()))
            # minV = np.min((minV,val.min()))

            valS.append(val)
            crdS.append(crd)

            iprint("Successfully processed : {}".format(pth))
        except KeyboardInterrupt:
            iprint("User exit")
            sys.exit(-1)
        # except:
        #     eprint("Failed to process : {}".format(pth))

    iprint("Data read : COMPLETE")
    for pth,crd,val in zip(args.count_files,crdS,valS):
        if args.threshold is not None:
            iprint("using cutoff : {}".format(args.threshold))
            val[val < args.threshold] = 0.0
            _threshold = args.threshold
        else:
            _threshold = "None"


        edgecolor = np.zeros((val.shape[0],4))
        edgecolor[:,3] = 0.4

        figsize = (10,10)
        fig,ax = plt.subplots(1,1,
                            figsize = figsize)

        minV = eval(str(_threshold))

        ax.scatter(crd[:,0],
                crd[:,1],
                s = args.marker_size,
                c = val,
                vmin = minV,
                vmax = (maxV if maxV > minV else 1),
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

if __name__ == "__main__":
    main()
