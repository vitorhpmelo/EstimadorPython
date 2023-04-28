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
# conv=load_flow_FACTS_2(graph,prt=1)

#%%
conv=load_flow_FACTS(graph,inici=-1,prt=1,itmax=40)

#%%
for key,r in ramTCSC.items():
    print("{:e}".format(r.xtcsc))

ram.update(ramTCSC)
save_DMED_fp(graph,ram,sys)

#%%
state_ref=get_state(graph)

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)


#%%

print("FACTS1")

with open("conds.csv","w") as f:
    f.write("Estimador 1 \n")
SS_WLS_FACTS(graph,dfDMEDsr,ind_i,flatstart=4,pirntits=1,printcond=1)
#%%
print("FACTS2")
with open("conds.csv","a") as f:
    f.write("Estimador 2 \n")
SS_WLS_FACTS_2(graph,dfDMEDsr,ind_i,flatstart=5,pirntits=1,printcond=1)

#%%
#%%
# sys="IEEE14_alt"

# dfDBAR,dfDBRAN,dfDMED,dfDFACTS = read_files(sys)


# [bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)

# [ram,nbran]=create_bran(dfDBRAN,ind_i)

# graph=create_graph(bars,ram)
# conv = load_flow(graph,tol=1e-7) 

# save_DMED_fp(graph,ram,sys)



# #%%
# #%% fluxo de potência
# conv = load_flow(graph,tol=1e-7)#parei aqui, pensar se a Jacobiana do fluxo vai calcular as derivadas em relação a fluxo de potência e concatenar as matrizes, depois testar
# #%% salva o DMEDFP com todas as grandezas
# save_DMED_fp(graph,ram,sys)

# #%% Estimador com QR

# SS_WLS(graph,dfDMED,ind_i,solver="QR")
# #%% Estimador com Normal
# SS_WLS(graph,dfDMED,ind_i,solver="Normal")
# #%% Estimador Lagrangeano
# SS_WLS_lagrangian(graph,dfDMED,ind_i)


# # %%
