#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
from classes import *
from readfiles import *
from networkstruc import *
from SS import *
from meas_sampl import *
import pandas as pd
import numpy as np
from networkcalc import *
from BadData import *
import numpy.linalg as liang
import scipy.sparse.linalg as sliang 


#%% Lê arquivos e constroi a estrutura da rede

sys="IEEE14_tcsc"

dfDBAR,dfDBRAN,dfDMED,dfDFACTS=read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)
#%%
[ramTCSC,nbranTCSC]=create_TCSC(dfDFACTS,ind_i)
#%%


graph=create_graph(bars,ram)
#%%
addFACTSingraph(graph,ramTCSC)

#%%
Vinici_lf(graph)
zPf,var_x = create_z_x_loadflow_TCSC(graph)
[z,var_t,var_v]=create_z_x_loadflow(graph)
z=z+zPf

# implementação TCSC

#%%
Vinici_lf(graph)
dz=np.zeros(len(z))
H=np.zeros((len(z),len(var_t)+len(var_v)))
Hx=np.zeros((len(z),len(var_x)))
it=0
while it<20:
    calc_dz(z,graph,dz)
    calc_H_fp(z,var_t,var_v,graph,H)
    calc_H_fp_TCSC(z,var_x,graph,Hx)
    HTCSC=np.concatenate((H,Hx),axis=1)
    A=sparse.csc_matrix(HTCSC, dtype=float)
    dx=sliang.spsolve(A,dz)
    new_X(graph,var_t,var_v,dx)
    new_X_TCSCC(graph,len(var_t)+len(var_v),var_x,dx)
    if np.max(np.abs(dx))< 1e-8 and np.max(np.abs(dz)) <1e-10:
        print(it)
        prt_state(graph)
        break
    it=it+1
#%%
ram.update(ramTCSC)

save_DMED_fp(graph,ram,sys)

#%%




sys="IEEE14_alt"

dfDBAR,dfDBRAN,dfDMED,dfDFACTS = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)

[ram,nbran]=create_bran(dfDBRAN,ind_i)

graph=create_graph(bars,ram)
conv = load_flow(graph,tol=1e-7) 

save_DMED_fp(graph,ram,sys)



#%%
#%% fluxo de potência
conv = load_flow(graph,tol=1e-7)#parei aqui, pensar se a Jacobiana do fluxo vai calcular as derivadas em relação a fluxo de potência e concatenar as matrizes, depois testar
#%% salva o DMEDFP com todas as grandezas
save_DMED_fp(graph,ram,sys)

#%% Estimador com QR

SS_WLS(graph,dfDMED,ind_i,solver="QR")
#%% Estimador com Normal
SS_WLS(graph,dfDMED,ind_i,solver="Normal")
#%% Estimador Lagrangeano
SS_WLS_lagrangian(graph,dfDMED,ind_i)


# %%
