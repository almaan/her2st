#!/usr/bin/env python3

from gprofiler import GProfiler
import pandas as pd

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


res_pth = "/tmp/tls/tls-associated-top.tsv"
res = pd.read_csv(res_pth, sep = '\t',header = 0, index_col = 0)

organism = "hsapiens"


genes = res.index.values.tolist()

gp = GProfiler(return_dataframe = True)
res = gp.profile(organism = organism,
                 query = genes,
                 )

gobps = res.loc[res['source'].values == "GO:BP",:]
# is_sig = gobps['significant'].values
gobps = gobps.loc[gobps['significant'].values,:]

gobps.to_csv("/tmp/tls-enrich.tsv",sep = '\t',header = True,index =True)

label_dict = {"fontsize":20,"family":"calibri"}
standard_dict = {'family':'calibri'}
cmap = plt.cm.YlGnBu
n_top = 20
fig,ax = plt.subplots(1,1, figsize = (9,8))
y_pos = np.arange(n_top)[::-1] + 0.5
p_vals = gobps['p_value'].values[0:n_top]
x_pos = -np.log10(p_vals)

cvals = stretch(x_pos,0,cmap.N) / 256
cvals = cmap(cvals)

delta = -np.diff(x_pos)[0] 
ax.barh(y_pos,x_pos,color = cvals,edgecolor = 'black')

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
fig.savefig("/tmp/tls-enrich.png")



#print(gobps['name'].values.tolist())
