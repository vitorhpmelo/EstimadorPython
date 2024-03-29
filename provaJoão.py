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


#%%


sys="IEEE14"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)

#%%

SS_WLS(graph,dfDMED,ind_i,solver="Normal",tol=1e-5,prec_virtual=1e-5)
Cov=calcCovRes(graph,dfDMED,ind_i)
dfRe=renorm(graph,dfDMED,ind_i,Cov)
print("Maior Residuo")
print(dfRe[dfRe.Rn==max(dfRe.Rn)])

dfRe.to_csv("Residuos.csv")
#%%