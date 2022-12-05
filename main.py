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


#%% Constroi lê a estrutura da rede

prec_virtual=1e-7
sys="IEEE118"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)

[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)


#%% fluxo de potência
conv = load_flow(graph,tol=1e-7)
#%% salva o DMEDFP
save_DMED_fp(graph,ram,sys)

#%% Estimador com QR
#Rodar o EE
SS_WLS(graph,dfDMED,ind_i,solver="QR")


#%% Estimador gradientes conjugados

SS_WLS(graph,dfDMED,ind_i,solver="cg")

# %% Estimador Lagrangeano

SS_WLS_lagrangian(graph,dfDMED,ind_i)
# %% Residuos normalizados
Cov=calcCovRes(graph,dfDMED,ind_i)
dfRe=renorm(graph,dfDMED,ind_i,Cov)
dfRe.to_csv("Residuos.csv")



#%%

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED("IEEE14",prec,graph,ram)
dfDMEDr=insert_res(dfDMEDsr)
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR")
SS_WLS_lagrangian(graph,dfDMEDsr,ind_i)
Cov=calcCovRes(graph,dfDMEDr,ind_i)
dfRe=renorm(graph,dfDMEDr,ind_i,Cov)
dfRe.to_csv("Residuos.csv")
# %%
