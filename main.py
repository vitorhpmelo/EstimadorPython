#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%

from classes import *
from readfiles import *
from networkstruc import *
import pandas as pd
import numpy as np
from networkcalc import *
import scipy.sparse.linalg as sliang 
import scipy.sparse as sparse 

#%%


sys="IEEE14"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)

[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)




# %%

# Vinici(graph)
Vinici_lf(graph)

# %%
# PowerFlows(ram,graph,print=1)
# PowerInjc(graph,print=1)

varlist=[]
# %%

[z,var_t,var_v]=create_z_x_loadflow(graph)
# %%


# %%
dx=np.ones(len(var_t))
it=0
while(np.amax(np.abs(dx))>1e-9 and it <20):
    print(np.amax(np.abs(dx)))
    dz=calc_dz(z,graph)
    H=calc_H(z,var_t,var_v,graph)
    A=sparse.csc_matrix(H, dtype=float)
    dx=sliang.spsolve(A,dz)
    new_X(graph,var_t,var_v,dx)
    it=it+1


# %%

for no in graph:
    s="Barra: {:d} | V : {:f} | teta : {:f}".format(no.bar.id,no.V,no.teta*180/np.pi)
    print(s)
# %%
