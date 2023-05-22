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

sys="IEEE14"

#% 5% 10% 15% 20%


# sample=list(np.random.randint(low=1, high=20, size=4, dtype=int))
sample=[12, 2, 19, 17]
# sample=[184, 117, 67, 127, 84, 135, 55, 29, 159, 139, 90, 26, 36, 88, 114, 113, 91, 69, 193, 49, 66, 158, 84, 86, 104, 40, 145, 64, 90, 1, 31, 92, 152, 72, 127, 27, 17, 147]
# sample = [128,  26,  69, 167,  16,  65,  78,  91, 147,  32,  39, 170, 177, 103, 163,  19,  81,  21, 111, 123,  41,  63,  13,  35,  99,  46,161, 125, 132, 184,5,70,44, 28 ,85,97,64]
# sample = [124, 89, 37, 22, 175, 101, 31, 184, 187, 136, 17, 195, 10, 183, 100, 50, 154, 63, 60, 8, 130, 170, 119, 1, 152, 85, 139, 46, 189, 23, 177, 4, 190, 72, 193, 53, 105, 160]
#%%
# sample = [128,  26,  69, 167,  16,  65,  78,  91, 147,  32,  39, 170, 177, 103, 163,  19,  81,  21, 111, 123,  41,  63,  13,  35,  99,  46,161]
# sample = [128,  26,  69, 167,  16,  65,  78,  91, 147,  32,  39, 170, 177, 103, 163,  19,  81,  21, 111]
# sample = [128,  26,  69, 167,  16,  65,  78,  91, 147]

#%%

dfDBAR,dfDBRAN,dfDMED = read_files(sys)
dfDBRANref=dfDBRAN.copy()

#%%
with open("conds.csv","w") as file:
    file.write("Inicio")

with open("iter.csv","w") as file:
    file.write("Inicio")
     
        
fatlist=[0.80,0.40,0.20,0,-0.20,-0.40,-0.80]
lstIEEE14=[3,2,1,0]
lstIEEE118=[27,18,9,0]
for value in lstIEEE14:
    branches=sample[value:]
    for fat in fatlist:
        dfDBRAN=dfDBRANref.copy()
        for bra in branches:
            dfDBRAN.loc[dfDBRAN["id"]==bra,"r"]= (1-fat)*dfDBRAN.loc[dfDBRAN["id"]==bra,"r"]
            dfDBRAN.loc[dfDBRAN["id"]==bra,"x"]= (1-fat)*dfDBRAN.loc[dfDBRAN["id"]==bra,"x"]
        with open("conds.csv","a") as file:
            file.write("value{},fat{},".format(value,fat))
        with open("iter.csv","a") as file:
            file.write("\nvalue{},fat{},".format(value,fat))


        [bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
        [ram,nbran]=create_bran(dfDBRAN,ind_i)
        #%%

        network=netinfo(nbars,nbran,2*nbars-1,nteta=nbars-1,nv=nbars)

        graph=create_graph(bars,ram)


        conv = load_flow(graph,tol=1e-7) #possivel erro, inicialização da referência

        if conv ==0:
            print("valor de impedancia incompantivel flux não convergiu")
            quit()
        state_ref=get_state(graph)

        save_DMED_fp(graph,ram,sys)

        prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5}
        dfDMEDsr=create_DMED(sys,prec,graph,ram)

        #%%

        print("Metodo QR")
        SS_WLS(graph,dfDMEDsr,ind_i,solver="QR",printmat=1,printcond=1,prinnormgrad=1)
        print("Metodo Normal")
        SS_WLS(graph,dfDMEDsr,ind_i,solver="Normal",printmat=1,printcond=1,prinnormgrad=1)
        print("Metodo Lagrangeano")
        SS_WLS_lagrangian(graph,dfDMED,ind_i,printmat=1,printcond=1,printnormgrad=1)

# %%
print(sample)