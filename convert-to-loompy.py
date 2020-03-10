#!/usr/bin/env python3

import pandas as pd
import numpy as np

import loompy as lp


import sys
import os.path as osp
import re

cnt_pth = sys.argv[1]
mta_pth = sys.argv[2]
out_dir = sys.argv[3]

# test suite
#cnt_pth = "/home/alma/Documents/PhD/papers/HER2/data/bc-alex/test.tsv"
#mta_pth = "/home/alma/Documents/PhD/papers/HER2/data/bc-alex/CTP_HER2_xFIVE_metadata.csv"
#out_dir = "/tmp" 

# load meta data
mta = pd.read_csv(mta_pth,
                  sep = ',',
                  index_col = 0,
                  header = 0,
                  )


# Inspect meta file
print(mta.iloc[0:10,0:3])

# load count data
cnt = pd.read_csv(cnt_pth,
                  sep = ',',
                  index_col = 0,
                  header = 0,
                  )

# Inspect count file
print(cnt.iloc[0:10,0:3])


# get interescting elements and order them

inter = cnt.columns.intersection(mta.index)

mta = mta.loc[inter,:]
cnt = cnt.loc[:,inter]

ct_major = mta['celltype_major'].values
ct_subset = mta['celltype_subset'].values
barcode = mta.index.values
gene = cnt.index.values

col_attrs = dict(CellID = barcode,
                MajorCellType = ct_major,
                SubsetCellType = ct_subset,
                )

row_attrs = dict(Gene = gene )

basename = re.sub('\\.tsv|\\.tsv\\.gz|\\.csv','.loom',osp.basename(cnt_pth))
oname = osp.join(out_dir,basename)
ds = lp.create(oname,
               cnt.values,
               row_attrs,
               col_attrs)
