#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
from classes import *
from readfiles import *
from networkstruc import *
from SS import *
import pandas as pd
import numpy as np
from networkcalc import *
import numpy.linalg as liang
import scipy.sparse.linalg as sliang 





prec_virtual=1e-7
sys="IEEE14"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)

[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)
#%%
#conv = load_flow(graph,tol=1e-10)

#save_DMED_fp(graph,ram,sys)

#%%
SS_WLS(graph,dfDMED,ind_i)
# %%
[z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
W=create_W(z,flag_ones=0,prec_virtual=1e-5)
Wmei=np.sqrt(np.diag(W))

Vinici(graph,flatStart=1)
H=np.zeros((len(z),len(var_t)+len(var_v)))
dz=np.zeros(len(z))
calc_dz(z,graph,dz)
calc_H_EE(z,var_t,var_v,graph,H)


#R=NormalEQ(H,W,dz)

# %%
