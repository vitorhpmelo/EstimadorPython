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

sys="IEEE118_tcsc_2"
<<<<<<< HEAD
increase=20
flatstart=4
=======
increase=5
flatstart=5
>>>>>>> refs/remotes/origin/TCSC


dfDBAR,dfDBRAN,dfDMED,dfDFACTS=read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)
[ramTCSC,nbranTCSC]=create_TCSC(dfDFACTS,ind_i)
#%%


graph=create_graph(bars,ram)

#%%
addFACTSingraph(graph,ramTCSC)


# conv=load_flow_FACTS(graph,inici=-1,prt=1,itmax=40)
conv=load_flow_FACTS_2(graph,prt=1)

ram.update(ramTCSC)
save_DMED_fp(graph,ram,sys)

#%% Montagem do plano de medições

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)
dfDMEDr=insert_res(dfDMEDsr)
dfDMEDr.to_csv(sys+"/DMED.csv",header=None,index=None,float_format="%.9f")
#%%

#simulação de monte carlo

 #possivel erro, inicialização da referência
state_ref=get_state(graph)
state_FACTS_ref=get_state_TCSC(ramTCSC)

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)
N=100
nconvsA=[]
stateA=[]
stateA_X=[]
TemposTotaisA=[]
TemposiTA=[]
NumeroItsA=[]

nconvsB=[]
stateB_X=[]
stateB=[]
TemposTotaisB=[]
TemposiTB=[]
NumeroItsB=[]

pre1=1e-5
prec2=1e-8
i=0
<<<<<<< HEAD
itmax=10
=======
itmax=500
>>>>>>> refs/remotes/origin/TCSC
amostrasA=0
amostrasB=0
while((len(NumeroItsA)<101)or(len(NumeroItsB)<101)):
    dfDMEDr=insert_res(dfDMEDsr,i)
    print("{:d}/{:d}".format(i+1,N))
    if (len(NumeroItsA)<101):
        [conv,nIT,tits,tf]=SS_WLS_FACTS_clean(graph,dfDMEDr,ind_i,flatstart=flatstart,tol=1e-5,tol2=1e-4)
        nconvsA.append(conv)
        print(conv)
        amostrasA=amostrasA+1
        if conv == 1:
            stateA.append(get_state(graph))
            stateA_X.append(get_state_TCSC(ramTCSC))
            TemposTotaisA.append(tf)
            TemposiTA.append(np.mean(tits))
            NumeroItsA.append(nIT)
        del nIT,tits,tf
    if (len(NumeroItsB)<101):
        [conv,nIT,tits,tf]=SS_WLS_FACTS_2_clean(graph,dfDMEDr,ind_i,flatstart=flatstart,tol=1e-5,tol2=1e-4)
        nconvsB.append(conv)
        print(conv)
        amostrasB=amostrasB+1
        if conv == 1:
            stateB.append(get_state(graph))
            stateB_X.append(get_state_TCSC(ramTCSC))
            TemposTotaisB.append(tf)
            TemposiTB.append(np.mean(tits))
            NumeroItsB.append(nIT)
        del nIT,tits,tf
    i=i+1
    if i>itmax:
        break

<<<<<<< HEAD
=======

>>>>>>> refs/remotes/origin/TCSC
# %% calculo erros

erro_A_V=[]
erro_A_T=[]
for ste in stateA:
    erro_A_V.append(np.abs(state_ref.v-ste.v))
    erro_A_T.append(np.abs(state_ref.t-ste.t))

erro_A_X=[]
for ste_x in stateA_X:
    for key,x_ref in state_FACTS_ref.items():
        erro_A_X.append(np.abs(x_ref-ste_x[key]))

EMA_A_V=np.mean(np.mean(erro_A_V))
EMA_A_T=np.mean(np.mean(erro_A_T))
EMA_A_X=np.mean(np.mean(erro_A_X))
EMA_A_total=np.mean([EMA_A_V,EMA_A_T,EMA_A_X])    



erro_B_V=[]
erro_B_T=[]
for ste in stateB:
    erro_B_V.append(np.abs(state_ref.v-ste.v))
    erro_B_T.append(np.abs(state_ref.t-ste.t))



erro_B_X=[]
for ste_x in stateB_X:
    for key,x_ref in state_FACTS_ref.items():
        erro_B_X.append(np.abs(x_ref-ste_x[key]))


EMA_B_V=np.mean(np.mean(erro_B_V))
EMA_B_T=np.mean(np.mean(erro_B_T))
EMA_B_X=np.mean(np.mean(erro_B_X))
EMA_B_total=np.mean([EMA_B_V,EMA_B_T,EMA_B_X])    




try:
    maxiteB=max(NumeroItsB)
except:
    maxiteB=np.nan



try:
    maxiteA=max(NumeroItsA)
except:
    maxiteA=np.nan



#%% salva os resultados
dA={"Solver":"A","prec":increase,"EMA_V":EMA_A_V,"EMA_T":EMA_A_T,"EMA_X":EMA_A_X,"EMA_total":EMA_A_total,"Media Tempo total":np.mean(TemposTotaisA),"Media Tempo iteçoes":np.mean(TemposiTA),"Númeoro Médio de Iterações":np.mean(NumeroItsA),"Max itera":maxiteA,"Numero Convs":np.sum(nconvsA),"Namos":amostrasA}
dB={"Solver":"B","prec":increase,"EMA_V":np.float64(EMA_B_V),"EMA_T":np.float64(EMA_B_T),"EMA_X":EMA_B_X,"EMA_total":EMA_B_total,"Media Tempo total":np.mean(TemposTotaisB),"Media Tempo iteçoes":np.mean(TemposiTB),"Númeoro Médio de Iterações":np.mean(NumeroItsB),"Max itera":maxiteB,"Numero Convs":np.sum(nconvsB),"Namos":amostrasB}



#%% 
dfResults=pd.DataFrame([dA,dB])
dfResults.to_csv("restTCSC/res"+sys+"str_"+str(flatstart)+"incr_"+str(increase)+".csv")



# %%
