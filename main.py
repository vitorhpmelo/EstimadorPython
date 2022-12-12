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

conv = load_flow(graph,tol=1e-7) #possivel erro, inicialização da referência
state_ref=get_state(graph)

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)
N=100
sateNormal=[]
stateQR=[]
statelagran=[]
TemposTotaisNormal=[]
TemposTotaisQR=[]
TemposTotaisLagrange=[]
TemposiTNormal=[]
TemposiTQR=[]
TemposiTLagrange=[]
NumeroItsNormal=[]
NumeroItsQR=[]
NumeroItslagran=[]

for i in range(N):
    dfDMEDr=insert_res(dfDMEDsr,i)
    print("{:d}/{:d}".format(i+1,N))
    # print(dfDMEDr)
    [T,tits,conv,nits]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="Normal",prec_virtual=1e-5)
    if conv == 1:
        sateNormal.append(get_state(graph))
        TemposTotaisNormal.append(T)
        TemposiTNormal.append(np.mean(tits))
        NumeroItsNormal.append(nits)
    del nits,tits,T
    [T,tits,conv,nits]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="QR",prec_virtual=1e-5)
    if conv == 1:
        stateQR.append(get_state(graph))
        TemposTotaisQR.append(T)
        TemposiTQR.append(np.mean(tits))
        NumeroItsQR.append(nits)
    [T,tits,conv,nits]=SS_WLS_lagrangian_clean(graph,dfDMEDr,ind_i)
    if conv == 1:
        statelagran.append(get_state(graph))
        TemposTotaisLagrange.append(T)
        TemposiTLagrange.append(np.mean(tits))
        NumeroItslagran.append(nits)  
    del nits,tits,T
    
# %%

erroNormalV=[]
erroNormalT=[]
for ste in sateNormal:
    erroNormalV.append(np.abs(state_ref.v-ste.v))
    erroNormalT.append(np.abs(state_ref.t-ste.t))

EMA_Normal_V=np.mean(np.mean(erroNormalV))
EMA_Normal_T=np.mean(np.mean(erroNormalT))
EMA_Normal_total=np.mean([EMA_Normal_V,EMA_Normal_T])    

erroQRV=[]
erroQRT=[]
for ste in stateQR:
    erroQRV.append(np.abs(state_ref.v-ste.v))
    erroQRT.append(np.abs(state_ref.t-ste.t))    

EMA_QR_V=np.mean(np.mean(erroQRV))
EMA_QR_T=np.mean(np.mean(erroQRT))    
EMA_QR_total=np.mean([EMA_QR_V,EMA_QR_T])    


erroLagranV=[]
erroLagranT=[]
for ste in statelagran:
    erroLagranV.append(np.abs(state_ref.v-ste.v))
    erroLagranT.append(np.abs(state_ref.t-ste.t))

EMA_Lagran_V=np.mean(np.mean(erroLagranV))
EMA_Lagran_T=np.mean(np.mean(erroLagranT))   
EMA_Lagran_total=np.mean([EMA_Lagran_V,EMA_Lagran_T])     

#%%


#%%
dNormal={"Solver":"Normal","EMA_V":EMA_Normal_V,"EMA_T":EMA_Normal_T,"EMA_total":EMA_Normal_total,"Media Tempo total":np.mean(TemposTotaisNormal),"Media Tempo iteçoes":np.mean(TemposiTNormal),"Númeoro Médio de Iterações":np.mean(NumeroItsNormal)}
dQR={"Solver":"QR","EMA_V":EMA_QR_V,"EMA_T":EMA_QR_T,"EMA_total":EMA_QR_total,"Media Tempo total":np.mean(TemposTotaisQR),"Media Tempo iteçoes":np.mean(TemposiTQR),"Númeoro Médio de Iterações":np.mean(NumeroItsQR)}
dLagran={"Solver":"Lagran","EMA_V":EMA_Lagran_V,"EMA_T":EMA_Lagran_T,"EMA_total":EMA_Lagran_total,"Media Tempo total":np.mean(TemposiTLagrange),"Media Tempo iteçoes":np.mean(TemposiTLagrange),"Númeoro Médio de Iterações":np.mean(NumeroItslagran)}
#%%
dfResults=pd.DataFrame([dNormal,dQR,dLagran])



# %%
