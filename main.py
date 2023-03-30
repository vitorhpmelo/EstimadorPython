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

sys="IEEE14_tcsc"

dfDBAR,dfDBRAN,dfDMED,dfDFACTS=read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)
#%%
[ramTCSC,nbranTCSC]=create_TCSC(dfDFACTS,ind_i)
#%%


graph=create_graph(bars,ram)
#%%
addFACTSingraph(graph,ramTCSC)






#%%
# Vinici_lf(graph)
zPf,var_x = create_z_x_loadflow_TCSC(graph)
[z,var_t,var_v]=create_z_x_loadflow(graph)
z=z+zPf

#%%
zcc_conv=[]
zcc_facts=[]
for item in z:
    if item.type==0 or item.type==2:
        if item.type==2:
            k=item.k
            m=item.m
            if (str(k)+"-"+str(m)) in graph[k].bFACTS_adjk.keys() or (str(m)+"-"+str(k)) in graph[k].bFACTS_adjm.keys():
                zcc_facts.append(item)
            else:
                zcc_conv.append(item)
        else:
            zcc_conv.append(item)
#%%

Hcc=np.zeros((len(zcc_conv)+len(zcc_facts),len(graph)))

i=0
d_injzcc={}
for item in zcc_conv:
    if item.type==0:
        soma=0
        k=item.k
        d_injzcc[k]=i
        for key,ram in graph[k].adjk.items():
            if ram.type!=3:
                Hcc[i][ram.para]=-1/(ram.x)
                soma=soma+1/(ram.x)
        for key,ram in graph[k].adjm.items():
            if ram.type!=3:
                Hcc[i][ram.de]=-1/(ram.x)
                soma=soma+1/(ram.x)
        Hcc[i][k]=soma
    elif item.type==2:
        k=item.k
        m=item.m
        Hcc[i][k]=1/(ram.x)
        Hcc[i][m]=-1/(ram.x)
    i=i+1

Hccx=np.zeros((len(zcc_conv)+len(zcc_facts),len(var_x)))
for item in zcc_facts:
    k=item.k
    m=item.m
    Hcc[i][k]=1/item.val
    Hcc[i][m]=-1/item.val
    if str(k)+"-"+str(m) in var_x.keys():
        Hccx[i][var_x[str(k)+"-"+str(m)]]=-1
    elif str(m)+"-"+str(k) in var_x.keys():
        Hccx[i][var_x[str(k)+"-"+str(m)]]=-1
    i=i+1

b=np.zeros(len(zcc_conv)+len(zcc_facts))
i=0
for item in zcc_conv:
    b[i]=item.val
for item in zcc_facts:
    k=item.k
    m=item.m    
    if k in d_injzcc.keys():
        b[d_injzcc[k]]=b[d_injzcc[k]]-item.val
    elif m in d_injzcc.keys():
        b[d_injzcc[m]]=b[d_injzcc[m]]+item.val
A=np.concatenate((Hcc,Hccx),axis=1)

A_red=np.delete(A,0,1)

x=liang.solve(A_red,b)


exit()
# implementação TCSC

#%%
# conv=load_flow_FACTS_2(graph)

#%%
# conv=load_flow_FACTS(graph)

#%%
# save_DMED_fp(graph,ram,sys)


#%%
SS_WLS_FACTS(graph,dfDMED,ind_i,flatstart=-1)

#%%
SS_WLS_FACTS_2(graph,dfDMED,ind_i,flatstart=-1)

exit()
#%%
sys="IEEE14_alt"

dfDBAR,dfDBRAN,dfDMED,dfDFACTS = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)

[ram,nbran]=create_bran(dfDBRAN,ind_i)

graph=create_graph(bars,ram)
conv = load_flow(graph,tol=1e-7) 

save_DMED_fp(graph,ram,sys)



#%%
#%% fluxo de potência
conv = load_flow(graph,tol=1e-7)#parei aqui, pensar se a Jacobiana do fluxo vai calcular as derivadas em relação a fluxo de potência e concatenar as matrizes, depois testar
#%% salva o DMEDFP com todas as grandezas
save_DMED_fp(graph,ram,sys)

#%% Estimador com QR

SS_WLS(graph,dfDMED,ind_i,solver="QR")
#%% Estimador com Normal
SS_WLS(graph,dfDMED,ind_i,solver="Normal")
#%% Estimador Lagrangeano
SS_WLS_lagrangian(graph,dfDMED,ind_i)


# %%
