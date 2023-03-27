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

#% 5% 10% 15% 20%

branches=[1]


dfDBAR,dfDBRAN,dfDMED = read_files(sys)

#%%



fat=0.05
for bra in branches:
    dfDBRAN.loc[dfDBRAN["id"]==bra,"r"]= (1-fat)*dfDBRAN.loc[dfDBRAN["id"]==bra,"r"]
    dfDBRAN.loc[dfDBRAN["id"]==bra,"x"]= (1-fat)*dfDBRAN.loc[dfDBRAN["id"]==bra,"x"]






[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)
#%%

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)


conv = load_flow(graph,tol=1e-7) #possivel erro, inicialização da referência

if conv ==0:
    print("valor de impedancia incompantivel flux não convergiu")
    quit()
state_ref=get_state(graph)

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)



print("Metodo QR")
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",printmat=1,printcond=1)
print("Metodo Normal")
SS_WLS(graph,dfDMEDsr,ind_i,solver="Normal",printmat=1,printcond=1)
print("Metodo Lagrangeano")
SS_WLS_lagrangian(graph,dfDMED,ind_i,printmat=1,printcond=1)

# %%
