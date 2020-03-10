#!/usr/bin/env python3

import torch as t
from torch.utils.data import Dataset,DataLoader
from typing import NoReturn, Dict, ClassVar, List
import torch.nn as nn
from torch.autograd import Variable

import argparse as arp
import pandas as pd
import numpy as np
import copy

class SpotData(Dataset):

    def __init__(self,
                 cnt : t.tensor,
                 lbl : t.tensor,
                )-> NoReturn:

        self.cnt = cnt
        # Normalize to avoid batch effects
        self.cnt = self.cnt / t.sum(self.cnt, dim =  1).reshape(-1,1)
        self.lbl = lbl

        self.S = self.cnt.shape[0]
        self.G = self.cnt.shape[1]

    def __len__(self,) -> int:
        return self.S

    def __getitem__(self,
                    x : int,
                   ) -> Dict[t.tensor,t.tensor]:

        return { 'expr': self.cnt[x,:],
                 'label':self.lbl[x]}


class LogisticRegression(nn.Module):

    def __init__(self,
                 n_genes : int,
                 ) -> NoReturn:

        super(LogisticRegression,self).__init__()

        self.G = n_genes
        self.B = nn.Parameter(t.empty(self.G))
        self.B0 = nn.Parameter(t.empty(1))
        nn.init.normal_(self.B,0,1)
        nn.init.normal_(self.B0)

        self.loss = t.tensor(0)
        #self.bce = nn.functional.binary_cross_entropy

    def forward(self,
                x : t.tensor,
               ) -> t.tensor:


        theta = t.einsum('bg, g -> b',x,self.B) + self.B0
        y_pred = t.sigmoid(theta)

        return y_pred


def fit(dataset : ClassVar,
        n_epochs : int,
        batch_size : int,
       ) -> ClassVar:


    n_train = int(0.8 * dataset.S)
    n_val = dataset.S - n_train
    train_dataset, val_dataset = t.utils.data.random_split(dataset,
                                                           (n_train,n_val))

    train_batch_size = min((n_train,batch_size))
    val_batch_size = min((n_val,batch_size))


    train_dataloader = DataLoader(train_dataset,
                                  batch_size = train_batch_size,
                                  shuffle = True, # be careful
                                  )

    val_dataloader = DataLoader(val_dataset,
                                  batch_size = val_batch_size,
                                  shuffle = True, # be careful
                                  )


    model = LogisticRegression(dataset.G)
    optim = t.optim.Adam(model.parameters(), lr = 0.1)
    loss_fun = nn.functional.binary_cross_entropy
    best_val_acc = 0.0

    # TODO : Update to compute loss in fit
    # use accuracy rather than loss
    for epoch in range(n_epochs):
        total_loss = 0
        model.train()
        for idx, data in enumerate(train_dataloader):
            optim.zero_grad()
            x_value = Variable(data['expr'])
            y_true = data['label']
            y_pred = model(x_value)
            loss = loss_fun(y_pred.reshape(-1,),y_true.reshape(-1,))
            total_loss += loss.item()
            loss.backward()
            optim.step()

        print('\r',end = '')
        print(f" Epoch : {epoch} | Training Loss : {total_loss}",end = '')

        if epoch % 10 == 0:
            model.eval()
            val_acc = 0.0
            den = 0.0
            for idx,data in enumerate(val_dataloader):
                y_true = data['label']
                x_value = data['expr']
                y_pred = t.round(model(x_value))

                val_acc += t.sum(y_pred.reshape(-1,) == y_true.reshape(-1,))
                den += y_pred.shape[0]

            val_acc = float(val_acc) / float(den)
            if best_val_acc < val_acc:
                best_model = copy.deepcopy(model)

            print(f"\nEpoch : {epoch} | Validation Accuracy {val_acc}",)

    return best_model

def eval(dataset : ClassVar,
          model : ClassVar,
          batch_size : int,
         ) -> float :

    dataloader = DataLoader(dataset,
                             batch_size = batch_size,
                             shuffle = True, # be careful
                            )

    model.eval()
    test_acc = t.tensor(0.0)

    for idx,data in enumerate(dataloader):
        y_true = data['label']
        x_value = data['expr']
        y_pred = t.round(model(x_value))
        test_acc += t.sum((y_pred.reshape(-1,) == y_true.reshape(-1,)))


    test_acc = float(test_acc) / float(dataset.S)
    print(f"\nTest Accuracy {test_acc}")

    return test_acc

def  generate_synth(n_spots : int,
                   n_genes : int,
                  ):

     return t.round(t.rand((n_spots,n_genes)) * 10), t.round(t.rand(n_spots))

def read_file(pth : str,
             ) -> pd.DataFrame:

    tmp = pd.read_csv(pth,
                      sep = '\t',
                      header = 0,
                      index_col = 0,
                     )

    return tmp

def read_data(cnt_pths : List[str],
              lbl_pths : List[str],
              column : str = "tumor",
              pos_lab : str = "tumor",
             )-> Dict:

    # Sort to ensure match
    cnt_pths.sort()
    lbl_pths.sort()

    assert len(cnt_pths) == len(lbl_pths), \
            "Count and label files do not match"

    clist = []
    llist = []

    genes = []
    n_spots = 0

    for k,(cp,lp) in enumerate(zip(cnt_pths,lbl_pths)):
        print(f"Reading section {k+1} / {len(cnt_pths)}")
        print(f"count path : {cp}")
        print(f"label path : {lp}")

        ctmp = read_file(cp)
        ltmp = read_file(lp)
        ltmp = ltmp.loc[:,column]

        inter = ctmp.index.intersection(ltmp.index)
        ctmp = ctmp.loc[inter,:]
        ltmp = ltmp.loc[inter]


        onehot = np.array([1 if label == pos_lab else \
                  0 for label in ltmp.values])


        ltmp  = pd.DataFrame(onehot,
                             columns = ['label'],
                             index = inter)

        index = pd.Index(["s_" + str(k) + "_" + ival for \
                         ival in inter.values ])

        ctmp.index = index
        ltmp.index = index

        clist.append(ctmp)
        llist.append(ltmp)

        genes.append(set(ctmp.columns.values.tolist()))
        n_spots += ctmp.shape[0]

    genes = pd.Index(list(set.intersection(*genes)))


    joint_cnt = pd.DataFrame([],columns = genes)
    joint_lbl = pd.DataFrame([],columns = ["label"])

    for cp,lp in zip(clist,llist):
        joint_cnt = pd.concat([joint_cnt,cp.loc[:,genes]])
        joint_lbl = pd.concat([joint_lbl,lp])



    return {'counts' : joint_cnt,
            'labels' : joint_lbl}

def trial():
    cnt, lbl = generate_synth(100,50)
    ds = SpotData(cnt,lbl)
    fit(ds,n_epochs = 25,batch_size = 16)

def main():

    prs = arp.ArgumentParser()
    prs.add_argument('-trc','--train_counts',
                     nargs = '+',
                     required = True,
                     help = ' '.join(['training count'
                                    'matrices',]
                                   )
                    )

    prs.add_argument('-tsc','--test_counts',
                     nargs = '+',
                     required = True,
                     help = ' '.join(['test count'
                                      'matrices',]
                                   )
                    )

    prs.add_argument('-trl','--train_labels',
                     nargs = '+',
                     required = True,
                     help = ' '.join(['training labels'
                                      '',]
                                   )
                    )

    prs.add_argument('-tsl','--test_labels',
                     nargs = '+',
                     required = True,
                     help = ' '.join(['test labels'
                                      '',]
                                   )
                    )

    prs.add_argument('-e','--epochs',
                     default = 100,
                     type = int,
                     help = 'number of epochs'
                    )

    prs.add_argument('-bs','--batch_size',
                     default = 128,
                     type = int,
                     help = '',
                    )

    args = prs.parse_args()

    train_data = read_data(args.train_counts,
                           args.train_labels,
                           )


    test_data = read_data(args.test_counts,
                          args.test_labels,
                         )

    intergenes = train_data['counts'].columns.intersection(test_data['counts'].columns)

    train_data['counts'] = train_data['counts'].loc[:,intergenes]
    test_data['counts'] = test_data['counts'].loc[:,intergenes]



    train_dataset = SpotData(t.tensor(train_data['counts'].values.astype(np.float32)),
                             t.tensor(train_data['labels'].values.astype(np.float32)),
                            )

    test_dataset = SpotData(t.tensor(test_data['counts'].values.astype(np.float32)),
                            t.tensor(test_data['labels'].values.astype(np.float32)),
                            )

    model = fit(train_dataset,
                n_epochs = args.epochs,
                batch_size = args.batch_size,
                )

    acc = eval(test_dataset,
               model,
               batch_size = args.batch_size)



if __name__ == '__main__':
    main()


