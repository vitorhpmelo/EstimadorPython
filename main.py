#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
from classes import *
from readfiles import *
from networkstruc import *
import pandas as pd
import numpy as np
from networkcalc import *


prec_virtual=1e-7
sys="IEEE14"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)

[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)
#%%
conv = load_flow(graph,tol=1e-10,prt=1)

# save_DMED_fp(graph,ram,sys)

#%%
[z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)


W=create_W(z,flag_ones=1)
Vinici(graph,flatStart=1)
H=np.zeros((len(z),len(var_t)+len(var_v)))
dz=np.zeros(len(z))
#%%

#%%
it=0
tol=1e-6
while(it <5):
    calc_dz(z,graph,dz)
    calc_H(z,var_t,var_v,graph,H)
    if(it==0 and 1):
            np.savetxt("HEE.csv",H,delimiter=",")
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    A=sparse.csc_matrix(G, dtype=float)
    dx=sliang.spsolve(A,grad)
    new_X(graph,var_t,var_v,dx)
    print(np.amax(np.abs(dx)))
    if (np.amax(np.abs(dx))<tol):
        conv=1
        txt="Convergiu em {:d} iteracoes".format(it)
        print(txt)
        prt_state(graph)
        break

    it=it+1

# %%
