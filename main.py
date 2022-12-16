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

sys="IEEE14"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)


#%%
YBus=np.zeros([len(graph),len(graph)],dtype=complex)

for key,item in ram.items():
    i=item.de
    j=item.para
    YBus[i][j]=YBus[i][j]+item.Y[0][1]
    YBus[j][i]=YBus[j][i]+item.Y[1][0]
    YBus[i][i]=YBus[i][i]+item.Y[0][0]
    YBus[j][j]=YBus[j][j]+item.Y[1][1]

for no in graph:
    if no.FlagBS==1:
        i=no.id
        YBus[i][i]=YBus[i][i]+complex(0,no.Bs)

np.savetxt(sys+"Ybus.csv",YBus)


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
