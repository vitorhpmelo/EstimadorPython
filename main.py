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
sateNormal5=[]
sateNormal7=[]
stateQR5=[]
stateQR7=[]
statelagran=[]
TemposTotaisNormal5=[]
TemposTotaisNormal7=[]
TemposTotaisQR5=[]
TemposTotaisQR7=[]
TemposTotaisLagrange=[]
TemposiTNormal5=[]
TemposiTNormal7=[]
TemposiTQR5=[]
TemposiTQR7=[]
TemposiTLagrange=[]
NumeroItsNormal5=[]
NumeroItsNormal7=[]
NumeroItsQR5=[]
NumeroItsQR7=[]
NumeroItslagran=[]
nconvsNormal5=[]
nconvsNormal7=[]
nconvsQR5=[]
nconvsQR7=[]
nconvsLagran=[]
pre1=1e-5
prec2=1e-8
for i in range(N):
    dfDMEDr=insert_res(dfDMEDsr,i)
    print("{:d}/{:d}".format(i+1,N))
    # print(dfDMEDr)
    [T,tits,conv,nits]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="Normal",prec_virtual=pre1)
    nconvsNormal5.append(conv)
    if conv == 1:
        sateNormal5.append(get_state(graph))
        TemposTotaisNormal5.append(T)
        TemposiTNormal5.append(np.mean(tits))
        NumeroItsNormal5.append(nits)
    del nits,tits,T
    [T,tits,conv,nits]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="QR",prec_virtual=pre1)
    nconvsQR5.append(conv)
    if conv == 1:
        stateQR5.append(get_state(graph))
        TemposTotaisQR5.append(T)
        TemposiTQR5.append(np.mean(tits))
        NumeroItsQR5.append(nits)
    [T,tits,conv,nits]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="Normal",prec_virtual=prec2)
    nconvsNormal7.append(conv)
    if conv == 1:
        sateNormal7.append(get_state(graph))
        TemposTotaisNormal7.append(T)
        TemposiTNormal7.append(np.mean(tits))
        NumeroItsNormal7.append(nits)
    del nits,tits,T
    [T,tits,conv,nits]=SS_WLS_clean(graph,dfDMEDr,ind_i,solver="QR",prec_virtual=prec2)
    nconvsQR7.append(conv)
    if conv == 1:
        stateQR7.append(get_state(graph))
        TemposTotaisQR7.append(T)
        TemposiTQR7.append(np.mean(tits))
        NumeroItsQR7.append(nits)    
    [T,tits,conv,nits]=SS_WLS_lagrangian_clean(graph,dfDMEDr,ind_i)
    nconvsLagran.append(conv)
    if conv == 1:
        statelagran.append(get_state(graph))
        TemposTotaisLagrange.append(T)
        TemposiTLagrange.append(np.mean(tits))
        NumeroItslagran.append(nits)  
    del nits,tits,T
    
# %% calculo erros

erroNormal5V=[]
erroNormal5T=[]
for ste in sateNormal5:
    erroNormal5V.append(np.abs(state_ref.v-ste.v))
    erroNormal5T.append(np.abs(state_ref.t-ste.t))

EMA_Normal_5_V=np.mean(np.mean(erroNormal5V))
EMA_Normal_5_T=np.mean(np.mean(erroNormal5T))
EMA_Normal_5_total=np.mean([EMA_Normal_5_V,EMA_Normal_5_T])    

erroQR5V=[]
erroQR5T=[]
for ste in stateQR5:
    erroQR5V.append(np.abs(state_ref.v-ste.v))
    erroQR5T.append(np.abs(state_ref.t-ste.t))    

EMA_QR_5_V=np.mean(np.mean(erroQR5V))
EMA_QR_5_T=np.mean(np.mean(erroQR5T))    
EMA_QR_5_total=np.mean([EMA_QR_5_V,EMA_QR_5_T])    


erroNormal7V=[]
erroNormal7T=[]
for ste in sateNormal7:
    erroNormal7V.append(np.abs(state_ref.v-ste.v))
    erroNormal7T.append(np.abs(state_ref.t-ste.t))

EMA_Normal_7_V=np.mean(np.mean(erroNormal7V))
EMA_Normal_7_T=np.mean(np.mean(erroNormal7T))
EMA_Normal_7_total=np.mean([EMA_Normal_7_V,EMA_Normal_7_T])    

erroQR7V=[]
erroQR7T=[]
for ste in stateQR7:
    erroQR7V.append(np.abs(state_ref.v-ste.v))
    erroQR7T.append(np.abs(state_ref.t-ste.t))    

EMA_QR_7_V=np.mean(np.mean(erroQR7V))
EMA_QR_7_T=np.mean(np.mean(erroQR7T))    
EMA_QR_7_total=np.mean([EMA_QR_7_V,EMA_QR_7_T])    



erroLagranV=[]
erroLagranT=[]
for ste in statelagran:
    erroLagranV.append(np.abs(state_ref.v-ste.v))
    erroLagranT.append(np.abs(state_ref.t-ste.t))

EMA_Lagran_V=np.mean(np.mean(erroLagranV))
EMA_Lagran_T=np.mean(np.mean(erroLagranT))   
EMA_Lagran_total=np.mean([EMA_Lagran_V,EMA_Lagran_T])     



#%%
dNormal5={"Solver":"Normal","prec":pre1,"EMA_V":EMA_Normal_5_V,"EMA_T":EMA_Normal_5_T,"EMA_total":EMA_Normal_5_total,"Media Tempo total":np.mean(TemposTotaisNormal5),"Media Tempo iteçoes":np.mean(TemposiTNormal5),"Númeoro Médio de Iterações":np.mean(NumeroItsNormal5),"Max itera":max(NumeroItsNormal5),"Numero Convs":np.sum(nconvsNormal5)}
dQR5={"Solver":"QR","prec":pre1,"EMA_V":EMA_QR_5_V,"EMA_T":EMA_QR_5_T,"EMA_total":EMA_QR_5_total,"Media Tempo total":np.mean(TemposTotaisQR5),"Media Tempo iteçoes":np.mean(TemposiTQR5),"Númeoro Médio de Iterações":np.mean(NumeroItsQR5),"Max itera":max(NumeroItsQR5),"Numero Convs":np.sum(nconvsQR5)}
dNormal7={"Solver":"Normal","prec":prec2,"EMA_V":EMA_Normal_7_V,"EMA_T":EMA_Normal_7_T,"EMA_total":EMA_Normal_7_total,"Media Tempo total":np.mean(TemposTotaisNormal7),"Media Tempo iteçoes":np.mean(TemposiTNormal7),"Númeoro Médio de Iterações":np.mean(NumeroItsNormal7),"Max itera":max(NumeroItsNormal7),"Numero Convs":np.sum(nconvsNormal7)}
dQR7={"Solver":"QR","prec":prec2,"EMA_V":EMA_QR_7_V,"EMA_T":EMA_QR_7_T,"EMA_total":EMA_QR_7_total,"Media Tempo total":np.mean(TemposTotaisQR7),"Media Tempo iteçoes":np.mean(TemposiTQR7),"Númeoro Médio de Iterações":np.mean(NumeroItsQR7),"Max itera":max(NumeroItsQR7),"Numero Convs":np.sum(nconvsQR7)}
dLagran={"Solver":"Lagran","prec":0,"EMA_V":EMA_Lagran_V,"EMA_T":EMA_Lagran_T,"EMA_total":EMA_Lagran_total,"Media Tempo total":np.mean(TemposTotaisLagrange),"Media Tempo iteçoes":np.mean(TemposiTLagrange),"Númeoro Médio de Iterações":np.mean(NumeroItslagran),"Max itera":max(NumeroItslagran),"Numero Convs":np.sum(nconvsLagran)}
#%%
dfResults=pd.DataFrame([dNormal5,dQR5,dNormal7,dQR7,dLagran])
dfResults.to_csv("resultados/"+sys+".csv")



# %%
