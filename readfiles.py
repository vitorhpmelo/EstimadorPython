from classes import *
import pandas as pd
import numpy as np

def read_files(sys):
    try:
        dfDBAR=pd.read_csv(sys+"/DBAR.csv",header=None,dtype={0:np.int64,1:np.int64})
        dfDBAR.columns=["id","type","V","teta","Pg","Qg","Pd","Qd","Bs"]
    except:
        print("Error while reading DBAR file")
        exit(1)

    try:
        dfDBRAN=pd.read_csv(sys+"/DBRAN.csv",header=None,dtype={0:np.int64,1:np.int64,2:np.int64,3:np.int64})
        dfDBRAN.columns=["id","type","de","para","r","x","bsh","tap"]
    except:
        print("Error while reading DBRAN file")
        exit(1)
    try:
        dfDMED=pd.read_csv(sys+"/DMED.csv",header=None)
        dfDMED.columns=["type","de","para","zmed","prec"]
    except:
        print("There is no DMED")
        dfDMED=[]
    return dfDBAR,dfDBRAN,dfDMED


def prt_state(graph):
    for no in graph:
        s="Barra: {:d} | V : {:f} | teta : {:f}".format(no.bar.id,no.V,no.teta*180/np.pi)
        print(s)


def save_DMED_fp(graph,ram,sys):
    medidas=[]
    Pinj=[]
    Qinj=[]
    Vmod=[]
    Pkm=[]
    Pmk=[]
    Qkm=[]
    Qmk=[]
    for no in graph:
        linha=[0,no.bar.id,-1,no.P(graph),1]
        Pinj.append(linha)
        linha=[1,no.bar.id,-1,no.Q(graph),1]
        Qinj.append(linha)
        linha=[4,no.bar.id,-1,no.V,1]
        Vmod.append(linha)

    for key,r in ram.items():
        linha=[2,graph[r.de].bar.id,graph[r.para].bar.id,r.Pf(graph,0),1.0]
        linha2=[3,graph[r.de].bar.id,graph[r.para].bar.id,r.Qf(graph,0),1.0]
        Pkm.append(linha)
        Qkm.append(linha2)
        linha=[2,graph[r.para].bar.id,graph[r.de].bar.id,r.Pf(graph,1),1.0]
        linha2=[3,graph[r.para].bar.id,graph[r.de].bar.id,r.Qf(graph,1),1.0]
        Pmk.append(linha)
        Qmk.append(linha2)

    medidas=Pinj+Qinj+Pkm+Qkm+Pmk+Qmk+Vmod
    dfDMED=pd.DataFrame(medidas,columns=["type","de","para","zmed","pre"])
    dfDMED.to_csv(sys+"/DMED_fp.csv",index=False,float_format="%.7f",header=False)
