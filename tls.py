#!/usr/bin/env python3

import pandas as pd
import numpy as np

import argparse as arp
import os.path as osp

import statsmodels.api as sm
from scipy.ndimage import gaussian_filter

from typing import List


def iprint(s : str) -> None:
    print("[INFO] : {}".format(s))

def read_file(f : str,
              sep : str = '\t',
              )->pd.DataFrame:

    return pd.read_csv(f,header = 0,index_col = 0, sep = sep)

def joiner(pths : List[str],
           verbose : bool = True,
           file_type : str = '',
           )->pd.DataFrame:

    out = pd.DataFrame([])
    for k,pth in enumerate(pths):

        if verbose:
            print("\r[INFO] : Reading{} file [{}/{}] : {}".format(" " + file_type,
                                                                  k + 1,
                                                                  len(pths),
                                                                  pth),
                  end ='')

        out = pd.concat((out,read_file(pth)),sort = False)
        new_index = [str(k) + '_' + x for x in out.index.values]
        out.index = new_index

    out[out.isna()] = 0.0
    print("")

    return out

def get_inflection_point(times : np.ndarray,
                         sigma : float = 10,
                         min_genes : int = 1,
                         n_fail : int = 100,
                         )-> int:

    f_times = gaussian_filter(times,sigma)
    f_d2 = gaussian_filter(np.gradient(np.gradient(f_times)),sigma)
    first = np.argmax(f_d2 > 0)
    f_d2[0:first] = 1
    ipoint = np.argmax(f_d2 <= 0)

    if ipoint < min_genes:
        ipoint = np.min((f_times,n_fail))
        print("[WARNING] : Did not find inflection point"
              " using {} features instead.".format(ipoint))
    else:
        print("[INFO] : Using {} top features".format(ipoint))

    return ipoint

prs = arp.ArgumentParser()
aa = prs.add_argument

aa("-c",
   "--count_files",
   nargs = '+')

aa("-p",
   "--proportion_files",
   nargs = '+')

aa("-o",
   "--out_dir",
   default = '/tmp')

aa("-z",
   "--threshold",
   default = False,
   action = 'store_true',
   # type = int,
   )

aa("-t",
   "--tag",
   default = None,
   )

aa("-i","--use_intercept",
   default = False,
   action = 'store_true',
   )


aa("-a","--alpha",
   default = "0",
   type = str,
   )

args = prs.parse_args()


if args.tag is not None:
    tag = '-' + args.tag
else:
    tag = ''


out_pth = osp.join(args.out_dir,
                      "tls-associated-$$" + tag + ".tsv")


args.alpha = eval(args.alpha)

prop_pths = args.proportion_files
cnt_pths = args.count_files

prop_pths.sort()
cnt_pths.sort()

prop = joiner(prop_pths,file_type = 'proportion')
cnt = joiner(cnt_pths, file_type = 'count')


libSize = cnt.values.sum(axis = 1,keepdims= True)
nObs = cnt.values.sum(axis = 0, keepdims = True)

keep_genes = nObs.flatten() > 10
keep_spots = libSize.flatten() > 0

# std = cnt.values.std(axis=1,keepdims = True)
# mu = cnt.values.mean(axis=1,keepdims = True)

cnt = pd.DataFrame(np.divide(cnt.values, libSize,where = libSize > 0),
                   index = cnt.index,
                   columns = cnt.columns,
                   )

cnt = cnt.iloc[keep_spots,keep_genes]

inter = cnt.index.intersection(prop.index)

prop = prop.loc[inter,:]
cnt = cnt.loc[inter,:]


if args.use_intercept:
    iprint(">>> including intercept")
    cnt['intercept'] = np.ones(cnt.shape[0])

n_spots = cnt.shape[0]

type_1 = 'B-cells'
type_2 = 'T-cells'

pos_1 = np.argmax(prop.columns == type_1)
pos_2 = np.argmax(prop.columns == type_2)


jprod = np.zeros(n_spots)

iprint("computing TLS-Score")
for s in range(n_spots):
    vec = prop.values[s,:].reshape(-1,1)
    prod = np.dot(vec,vec.T)
    nprod = prod / prod.sum()
    N = prod.shape[0]

    jprod[s] = nprod[pos_1,pos_2] * 2
    jprod[s] -= (nprod.sum() / (0.5*(N**2 + N)))


jprod = pd.DataFrame(jprod,
                     index = cnt.index,
                     columns = ['probability'])


iprint("computing coefficients")
if args.alpha > 0 :

    iprint(">>> using regularization with : alpha = {}".format(args.alpha))

mod = sm.OLS(jprod.values,cnt.values)
res = mod.fit_regularized(L1_wt = 0,
                          alpha = args.alpha,
                          )

coefs = res.params

coefs = pd.Series(coefs,
                  index = cnt.columns,
                   )

ordr = np.argsort(coefs.values)[::-1]
coefs = coefs.iloc[ordr]


complete_out_pth = out_pth.replace("$$","complete")
iprint("saving complete results to : {}".format(complete_out_pth))
coefs.to_csv(complete_out_pth,
             sep = '\t',
             header = True,
             index = True)


if args.threshold:
    bs =  np.arange(coefs.shape[0])
    if args.use_intercept:
        b0 = np.argmax(coefs.index == 'intercept')
        bs = bs[coefs.index != 'intercept']

    pos = get_inflection_point(coefs.values[bs])
    sel = bs[0:pos]

    if args.use_intercept:
        sel = np.append(sel,b0)

    coefs = coefs.iloc[sel]

    top_out_pth = out_pth.replace("$$","top")
    iprint("saving top results to : {}".format(top_out_pth))
    coefs.to_csv(top_out_pth,
                sep = '\t',
                header = True,
                index = True)

iprint("TLS-analysis completed")
