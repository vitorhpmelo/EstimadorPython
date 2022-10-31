from csv import DictReader
from classes import *
import numpy as np

def creat_bar(dfDBAR):  
    bars=[]
    i=0
    nvar=0
    pv=[]
    pq=[]
    ind_id={}
    for idx, row in dfDBAR.iterrows():
        item=bar(int(row["id"]),int(row["type"]),i)
        item.V=row.V
        item.teta=row.teta*np.pi/180
        item.Pg=row["Pg"]
        item.Qg=row["Qg"]
        item.Pd=row["Pd"]
        item.Qd=row["Qd"]
        item.Bs=row["Bs"]
        if int(row["type"])==1:
            pv.append(i)
        elif int(row["type"])==2:
            pq.append(i)
        bars.append(item)
        ind_id[int(row["id"])]=i
        i=i+1
    return bars,i,pv,pq,ind_id

def create_bran(dfDBRAN,ind_i):
    ram={}
    i=0
    for id, row in dfDBRAN.iterrows():
        key=str(ind_i[int(row["de"])])+"-"+str(ind_i[int(row["para"])])
        item=branch(int(row["id"]),ind_i[int(row["de"])],ind_i[int(row["para"])],int(row["type"]),i)
        item.x=row["x"]
        item.r=row["r"]
        item.bsh=complex(0,row["bsh"]/2)
        item.tap=row["tap"]
        item.cykm()
        item.twoPortCircuit()
        ram[key]=item    
        i+=1
    return ram,i

def create_graph(bars,ram):
    graph=[]
    for item in bars:
        node=node_graph(item.i,item)
        if abs(item.Bs)>1e-8:
            node.FlagBS=1
            node.Bs=item.Bs
        graph.append(node)
    for key,item in ram.items():
        k=int(key.split("-")[0])
        m=int(key.split("-")[1])
        graph[k].adjk.update({key:item})
        graph[m].adjm.update({key:item})
    
    return graph