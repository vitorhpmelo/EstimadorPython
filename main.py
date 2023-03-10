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

sys="IEEE5TCSC"

dfDBAR,dfDBRAN,dfDMED,dfDFACTS = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)
[ramTCSC,nbranTCSC]=create_TCSC(dfDFACTS,ind_i)



graph=create_graph(bars,ram)
addFACTSingraph(graph,ramTCSC)

#%%
zPf,var_x = create_z_x_loadflow_TCSC(graph)

# implementação TCSC
#%%
dz=[0]
calc_dz(zPf,graph,dz)

#%%
Ybusmatlab=np.loadtxt("Ybusiee14.txt",delimiter=',',dtype=complex)
#%% fluxo de potência
conv = load_flow(graph,tol=1e-7)
#%% salva o DMEDFP com todas as grandezas
save_DMED_fp(graph,ram,sys)

#%% Estimador com QR

SS_WLS(graph,dfDMED,ind_i,solver="QR")
#%% Estimador com Normal
SS_WLS(graph,dfDMED,ind_i,solver="Normal")
#%% Estimador Lagrangeano
SS_WLS_lagrangian(graph,dfDMED,ind_i)


# %%
