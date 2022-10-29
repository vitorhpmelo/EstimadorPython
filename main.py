#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%

from classes import *
from readfiles import *
from networkstruc import *
import pandas as pd
import numpy as np




#%%


sys="IEEE14"

dfDBAR,dfDBRAN,dfDMED = read_files(sys)

[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

graph=create_graph(bars,ram)


# %%
