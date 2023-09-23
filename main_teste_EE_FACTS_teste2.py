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

sys="IEEE118_rakp2009"
per=60


dfDBAR,dfDBRAN,dfDMED,dfDFACTS=read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)
#%%
[ramTCSC,nbranTCSC]=create_TCSC(dfDFACTS,ind_i)

[busSVC,BUS_SVC]=create_SVC(dfDFACTS,ind_i)

[ramUPFC,nbranUPFC]=create_UPFC(dfDFACTS,ind_i)
#%%
for key,tcsc in ramTCSC.items():
    ramTCSC[key].xtcsc_ini=ramTCSC[key].xtcsc_ini*(1+per/100)

for key,_ in busSVC.items():
    busSVC[key].Bini=busSVC[key].Bini*(1+per/100)


for key,_ in ramUPFC.items():
    ramUPFC[key].Vse_ini=ramUPFC[key].Vse_ini*(1+per/100)
    ramUPFC[key].Vsh_ini=ramUPFC[key].Vsh_ini*(1+per/100)
    ramUPFC[key].t_se_ini=ramUPFC[key].t_se_ini*(1+per/100)
    ramUPFC[key].t_sh_ini=ramUPFC[key].t_sh_ini*(1+per/100)




#%%
graph=create_graph(bars,ram)

addTCSCingraph(graph,ramTCSC)

addSVCingraph(graph,busSVC)

addUPFCingraph(graph,ramUPFC)

#%%

print("Estimador 1")
print("FACTS with BC")
it1=SS_WLS_FACTS_withBC(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)
#%%
print("Estimador 2")
print("FACTS no BC")
it2=SS_WLS_FACTS_noBC(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)
#%%
print("Estimador 3")
print("FACTS with BC and jump first iterations facts")
it3=SS_WLS_FACTS_withBC_limalphavarfacts(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)

# %%
print("Estimador 4")
print("FACTS with BC lim for X and jump first iterations facts")
it4=SS_WLS_FACTS_withBC_limvarfacts(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)
#%%

print("Estimador 5")
print("FACTS with BC lim for X and jump first iterations facts")
it5=SS_WLS_FACTS_withBC_itvarfacts(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)
# %%

print("{}\n{}\n{}\n{}\n{}".format(it1,it2,it3,it4,it5))
# %%