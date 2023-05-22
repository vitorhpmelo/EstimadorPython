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


# sample = [128,  26,  69, 167,  16,  65,  78,  91, 147,  32,  39, 170, 177, 103, 163,  19,  81,  21, 111, 123,  41,  63,  13,  35,  99,  46,161, 125, 132, 184]
sample = [124, 89, 37, 22, 175, 101, 31, 184, 187, 136, 17, 195, 10, 183, 100, 50, 154, 63, 60, 8, 130, 170, 119, 1, 152, 85, 139, 46, 189, 23, 177, 4, 190, 72, 193, 53, 105, 160]

#%%


lstIEEE118=[27,18,9,0]
branches=sample
fat=0.80
qnt=0
for bra in branches[qnt:]:
    dfDBRAN.loc[dfDBRAN["id"]==bra,"r"]= (1-fat)*dfDBRAN.loc[dfDBRAN["id"]==bra,"r"]
    dfDBRAN.loc[dfDBRAN["id"]==bra,"x"]= (1-fat)*dfDBRAN.loc[dfDBRAN["id"]==bra,"x"]




#%% fluxo de potência
conv = load_flow(graph,tol=1e-7)
#%% salva o DMEDFP com todas as grandezas
save_DMED_fp(graph,ram,sys)



#%% Montagem do plano de medições

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)
dfDMEDr=insert_res(dfDMEDsr)
dfDMEDr.to_csv(sys+"/DMED.csv",header=None,index=None,float_format="%.9f")


#simulação de monte carlo

conv = load_flow(graph,tol=1e-7) #possivel erro, inicialização da referência
state_ref=get_state(graph)

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)
N=100
sateNormal5=[]
stateQR5=[]
statelagran=[]
TemposTotaisNormal5=[]
TemposTotaisQR5=[]
TemposTotaisLagrange=[]
TemposiTNormal5=[]
TemposiTQR5=[]
TemposiTLagrange=[]
NumeroItsNormal5=[]
NumeroItsQR5=[]
NumeroItslagran=[]
nconvsNormal5=[]
nconvsQR5=[]
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


erroLagranV=[]
erroLagranT=[]
for ste in statelagran:
    erroLagranV.append(np.abs(state_ref.v-ste.v))
    erroLagranT.append(np.abs(state_ref.t-ste.t))

EMA_Lagran_V=np.mean(np.mean(erroLagranV))
EMA_Lagran_T=np.mean(np.mean(erroLagranT))   
EMA_Lagran_total=np.mean([EMA_Lagran_V,EMA_Lagran_T])     



#%% salva os resultados
dNormal5={"Solver":"Normal","prec":fat,"EMA_V":EMA_Normal_5_V,"EMA_T":EMA_Normal_5_T,"EMA_total":EMA_Normal_5_total,"Media Tempo total":np.mean(TemposTotaisNormal5),"Media Tempo iteçoes":np.mean(TemposiTNormal5),"Númeoro Médio de Iterações":np.mean(NumeroItsNormal5),"Max itera":max(NumeroItsNormal5),"Numero Convs":np.sum(nconvsNormal5)}
dQR5={"Solver":"QR","prec":fat,"EMA_V":EMA_QR_5_V,"EMA_T":EMA_QR_5_T,"EMA_total":EMA_QR_5_total,"Media Tempo total":np.mean(TemposTotaisQR5),"Media Tempo iteçoes":np.mean(TemposiTQR5),"Númeoro Médio de Iterações":np.mean(NumeroItsQR5),"Max itera":max(NumeroItsQR5),"Numero Convs":np.sum(nconvsQR5)}
dLagran={"Solver":"Lagran","prec":fat,"EMA_V":EMA_Lagran_V,"EMA_T":EMA_Lagran_T,"EMA_total":EMA_Lagran_total,"Media Tempo total":np.mean(TemposTotaisLagrange),"Media Tempo iteçoes":np.mean(TemposiTLagrange),"Númeoro Médio de Iterações":np.mean(NumeroItslagran),"Max itera":max(NumeroItslagran),"Numero Convs":np.sum(nconvsLagran)}
#%% 
dfResults=pd.DataFrame([dNormal5,dQR5,dLagran])
dfResults.to_csv("resultados/"+sys+"_"+str(fat*100)+"_"+str(qnt)+".csv")



# %%
