#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
from classes import *
from readfiles import *
from networkstruc import *
import pandas as pd
import numpy as np
from networkcalc import *

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
[z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)


Vinici(graph,flatStart=1)
H=np.zeros((len(z),len(var_t)+len(var_v)))
dz=np.zeros(len(z))

#%%
W=create_W(z,flag_ones=0,prec_virtual=1e-5)
np.savetxt("W.csv",np.diag(W),delimiter=",")
#%%
it=0
tol=1e-5
while(it <10):
    calc_dz(z,graph,dz)
    calc_H_EE(z,var_t,var_v,graph,H)
    if(it==0):
        np.savetxt("H.csv",H,delimiter=",")
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    if(it==0):
        np.savetxt("G.csv",G,delimiter=",")
    A=sparse.csc_matrix(G)
    dx=sliang.spsolve(A,grad)
    #dx=np.linalg.solve(G,grad)
    new_X(graph,var_t,var_v,dx)
    print("max dx {:e} Cond G {:e}".format(np.amax(np.abs(dx)),np.linalg.cond(G)))
    if (np.amax(np.abs(dx))<tol):
        conv=1
        txt="Convergiu em {:d} iteracoes".format(it)
        print(txt)
        prt_state(graph)
        break
    it=it+1

# %%
