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

sys="IEEE14_rakp2009"


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


#%%
# print("Estimador 1")
# print("FACTS with BC")
# it3=SS_WLS_FACTS_withBC_itvarfacts(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)

# #%%
# print("Estimador 1")
# print("FACTS with BC")
it1=SS_WLS_FACTS_withBC(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)


#%%

it2=SS_WLS_FACTS_LM_BC(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)
# it2=SS_WLS_FACTS_grad(graph,dfDMED,ind_i,flatstart=2,pirntits=1,printcond=1,tol=1e-5,tol2=1e-4)

# %%
