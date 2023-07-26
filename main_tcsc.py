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

sys="5busAchaUPFC"


dfDBAR,dfDBRAN,dfDMED,dfDFACTS=read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)
#%%
[ramTCSC,nbranTCSC]=create_TCSC(dfDFACTS,ind_i)

[busSVC,BUS_SVC]=create_SVC(dfDFACTS,ind_i)

[ramUPFC,nbranUPFC]=create_UPFC(dfDFACTS,ind_i)




graph=create_graph(bars,ram)

addTCSCingraph(graph,ramTCSC)

addSVCingraph(graph,busSVC)

addUPFCingraph(graph,ramUPFC)

#%% confere derivadas Acha

# Vinici_DBAR(graph)
# ramUPFC["2-5"].dPpsdtp(graph)
# ramUPFC["2-5"].dQpsdtp(graph)
# ramUPFC["2-5"].dPpsdts(graph)
# ramUPFC["2-5"].dQpsdts(graph)
# ramUPFC["2-5"].dPpsdtse(graph)
# ramUPFC["2-5"].dQpsdtse(graph)
# ramUPFC["2-5"].dPpsdtsh(graph)
# ramUPFC["2-5"].dQpsdtsh(graph)
# ramUPFC["2-5"].dPspdtp(graph)
# ramUPFC["2-5"].dQspdtp(graph)
# ramUPFC["2-5"].dPspdts(graph)
# ramUPFC["2-5"].dQspdts(graph)
# ramUPFC["2-5"].dPspdtse(graph)
# ramUPFC["2-5"].dQspdtse(graph)
# ramUPFC["2-5"].dPspdtsh(graph)
# ramUPFC["2-5"].dQspdtsh(graph)

# #%%
# ramUPFC["2-5"].dPpsdVp(graph)
# ramUPFC["2-5"].dQpsdVp(graph)
# ramUPFC["2-5"].dPpsdVs(graph)
# ramUPFC["2-5"].dQpsdVs(graph)
# ramUPFC["2-5"].dPpsdVse(graph)
# ramUPFC["2-5"].dQpsdVse(graph)
# ramUPFC["2-5"].dPpsdVsh(graph)
# ramUPFC["2-5"].dQpsdVsh(graph)
# ramUPFC["2-5"].dPspdVp(graph)
# ramUPFC["2-5"].dQspdVp(graph)
# ramUPFC["2-5"].dPspdVs(graph)
# ramUPFC["2-5"].dQspdVs(graph)
# ramUPFC["2-5"].dPspdVse(graph)
# ramUPFC["2-5"].dQspdVse(graph)
# ramUPFC["2-5"].dPspdVsh(graph)
# ramUPFC["2-5"].dQspdVsh(graph)


# ramUPFC["2-5"].dIgdVp(graph)
# ramUPFC["2-5"].dIgdVs(graph)
# ramUPFC["2-5"].dIgdVse(graph)
# ramUPFC["2-5"].dIgdVsh(graph)

# ramUPFC["2-5"].dIgdtp(graph)
# ramUPFC["2-5"].dIgdts(graph)
# ramUPFC["2-5"].dIgdtse(graph)
# ramUPFC["2-5"].dIgdtsh(graph)



#%%
conv=load_flow_FACTS(graph,inici=1,prt=1,itmax=40)
#%%

try:
    for key,r in ramTCSC.items():
        print("{:e}".format(r.xtcsc))
except:
    pass

ram.update(ramTCSC)

save_DMED_fp(graph,ram,sys)



#%%
state_ref=get_state(graph)

prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
dfDMEDsr=create_DMED(sys,prec,graph,ram)


#%%
print("FACTS1")
dfDMEDr=insert_res(dfDMEDsr)
with open("conds.csv","w") as f:
    f.write("Estimador 1 \n")
SS_WLS_FACTS(graph,dfDMEDr,ind_i,flatstart=4,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)
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
