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

sys="IEEE118"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)




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
