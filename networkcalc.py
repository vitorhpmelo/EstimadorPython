from classes import *
import numpy as np
import pandas as pd

def Vinici(graph):
    '''
    Function to initate the voltages acording with DBAR
    '''
    for no in graph:
        no.V=no.bar.V
        no.teta=no.bar.teta


def PowerFlows(ram,graph,print=0):
    dpf={}
    dqf={}
    for key,bran in ram.items():
        k=key.split("-")[0]
        m=key.split("-")[1]
        dpf[k+"-"+m]=bran.Pf(graph,0)
        dpf[m+"-"+k]=bran.Pf(graph,1)
        dqf[k+"-"+m]=bran.Qf(graph,0)
        dqf[m+"-"+k]=bran.Qf(graph,1)
    if print == 1:
        dfFlows=pd.DataFrame(list(dpf.items()),columns=["from/to","PFlow"])
        dfFlows["QFlow"]=dfFlows["from/to"].map(dqf)
        dfFlows.to_csv("flows.csv",index=None)    
    return dpf,dqf

def PowerInjc(graph,print=0):
    dP={}
    dQ={}
    for item in graph:
        dP[item.id]=item.P(graph)
        dQ[item.id]=item.Q(graph)
    if print == 1:
        dfInj=pd.DataFrame(list(dP.items()),columns=["id","P"])
        dfInj["Q"]=dfInj["id"].map(dQ)    
        dfInj.to_csv("Inje.csv",index=None)  
         
    return dP,dQ 