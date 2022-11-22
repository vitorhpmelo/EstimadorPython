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




#%% Constroi lê a estrutura da rede

prec_virtual=1e-7
sys="IEEE30"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)

[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)
#%%
#Rodar o fluxo de potência

conv = load_flow(graph,tol=1e-10)

save_DMED_fp(graph,ram,sys)

#%%
#Rodar o EE
SS_WLS(graph,dfDMED,ind_i,solver="QR")


#%%


SS_WLS(graph,dfDMED,ind_i,solver="Normal")

# %%

SS_WLS_lagrangian(graph,dfDMED,ind_i)
# %%
