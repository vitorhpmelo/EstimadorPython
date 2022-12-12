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

sys="IEEE14"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)
N=2
sateNormal=[]
stateQR=[]
statelagran=[]
TemposTotaisNormal=[]
TemposTotaisQR=[]
TemposTotaisLagrange=[]
TemposiTNormal=[]
TemposiTQR=[]
TemposiTLagrange=[]
for i in range(N):
    dfDMEDr=insert_res(dfDMEDsr,i)
    print("{:d}/{:d}".format(i,N))
    # print(dfDMEDr)
    [T,tits,conv]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="Normal",prec_virtual=1e-5)
    if conv == 1:
        sateNormal.append(get_state(graph))
        TemposTotaisNormal.append(T)
        TemposiTNormal.append(np.mean(tits))
    [T,tits,conv]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="QR",prec_virtual=1e-5)
    if conv == 1:
        stateQR.append(get_state(graph))
        TemposTotaisQR.append(T)
        TemposiTQR.append(np.mean(tits))
    [T,tits,conv]=SS_WLS_lagrangian_clean(graph,dfDMEDr,ind_i)
    if conv == 1:
        sateNormal.append(get_state(graph))
        TemposTotaisLagrange.append(T)
        TemposiTLagrange.append(np.mean(tits))        
    


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
dfDMEDsr=create_DMED(sys,prec,graph,ram)
dfDMEDr=insert_res(dfDMEDsr)
dfDMEDr.to_csv(sys+"/DMED.csv",header=None,index=None,float_format="%.9f")
#%%
print("teste 5-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",prec_virtual=1e-5,printcond=1,prinnormgrad=1)
#%%
print("teste 1------------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="Normal",prec_virtual=1e-5,printcond=1,prinnormgrad=1)
print("teste 2-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="Normal",prec_virtual=1e-6,printcond=1,prinnormgrad=1)
print("teste 3-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="Normal",prec_virtual=1e-7,printcond=1,prinnormgrad=1)
print("teste 4-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="Normal",prec_virtual=1e-8,printcond=1,prinnormgrad=1)
print("teste 5-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="Normal",prec_virtual=1e-9,printcond=1,prinnormgrad=1)
print("teste 1------------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",prec_virtual=1e-5,printcond=1,prinnormgrad=1)
print("teste 2-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",prec_virtual=1e-6,printcond=1,prinnormgrad=1)
print("teste 3-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",prec_virtual=1e-7,printcond=1,prinnormgrad=1)
print("teste 4-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",prec_virtual=1e-8,printcond=1,prinnormgrad=1)
print("teste 5-----------------")
SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",prec_virtual=1e-9,printcond=1,prinnormgrad=1)

# print("teste 1------------------")
# SS_WLS(graph,dfDMEDsr,ind_i,solver="cg",prec_virtual=1e-5,printcond=1,prinnormgrad=1)
# print("teste 2-----------------")
# SS_WLS(graph,dfDMEDsr,ind_i,solver="cg",prec_virtual=1e-6,printcond=1,prinnormgrad=1)
# print("teste 3-----------------")
# SS_WLS(graph,dfDMEDsr,ind_i,solver="cg",prec_virtual=1e-7,printcond=1,prinnormgrad=1)
# print("teste 4-----------------")
# SS_WLS(graph,dfDMEDsr,ind_i,solver="cg",prec_virtual=1e-8,printcond=1,prinnormgrad=1)
# print("teste 5-----------------")
# SS_WLS(graph,dfDMEDsr,ind_i,solver="cg",prec_virtual=1e-9,printcond=1,prinnormgrad=1)
#%%
print("teste 6---Lagrangeano-------")
SS_WLS_lagrangian(graph,dfDMEDsr,ind_i,tol=1e-7,tol2=1e-5,printcond=1,printnormgrad=1,printmat=1)

#SS_WLS_lagrangian(graph,dfDMEDsr,ind_i)
#Cov=calcCovRes(graph,dfDMEDr,ind_i)
#dfRe=renorm(graph,dfDMEDr,ind_i,Cov)
#dfRe.to_csv("Residuos.csv")
# %%

#simulação de monte carlo


# %%
