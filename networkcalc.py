from classes import *
import numpy as np
import pandas as pd
from readfiles import *
import scipy.sparse.linalg as sliang 
import scipy.sparse as sparse 

def Vinici(graph):
    '''
    Function to initate the voltages acording with DBAR
    '''
    for no in graph:
        no.V=no.bar.V
        no.teta=no.bar.teta

def Vinici_lf(graph):
    '''
    Function to initate the voltages acording with DBAR
    '''
    for no in graph:
        if no.bar.type == 1 or no.bar.type == 0:
            no.V=no.bar.V
            no.teta=0
        else:
            no.V=1
            no.teta=0



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

def create_z_x_loadflow(graph):
    zP=[]
    zQ=[]
    var_t={}
    var_v={}
    i=0
    j=0
    for item in graph:
        if item.bar.type==1 or item.bar.type==2:
            mes=meas(item.id,-1,0,item.bar.Pg-item.bar.Pd,1)
            zP.append(mes)
            var_t[item.id]=i
            i=i+1
        if item.bar.type==2:
            mes=meas(item.id,-1,1,item.bar.Qg-item.bar.Qd,1)
            zQ.append(mes)
            var_v[item.id]=j
            j=j+1
    return zP+zQ,var_t,var_v



def calc_H(z,var_t,var_v,graph):
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    i=0
    n_teta=len(var_t)
    for item in z:
        soma1=0
        soma2=0
        if item.type==0:
            for key,branch in graph[item.k].adjk.items():#calcula as derivadas relativa as barras de e a do próprio angulo
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dPfdt(graph,0,item.k)
                if  branch.para in var_t.keys():
                    H[i][var_t[branch.para]]=branch.dPfdt(graph,0,branch.para)
            for key,branch in graph[item.k].adjm.items():
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dPfdt(graph,1,item.k) 
                if  branch.de in var_t.keys():
                    H[i][var_t[branch.de]]=branch.dPfdt(graph,1,branch.de)
            if  graph[item.k].bar.type!=0:
                H[i][var_t[item.k]]=soma1
            soma1=0
            for key,branch in graph[item.k].adjk.items(): ##derivadas modulo de tensão
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dPfdV(graph,0,item.k)
                if  branch.para in var_v.keys():
                    H[i][var_v[branch.para]+n_teta]=branch.dPfdV(graph,0,branch.para)
            for key,branch in graph[item.k].adjm.items():
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dPfdV(graph,1,item.k) 
                if  branch.de in var_v.keys():
                    H[i][var_v[branch.de]+n_teta]=branch.dPfdV(graph,1,branch.de)
            if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                H[i][var_v[item.k]+n_teta]=soma2
            soma2=0
        elif item.type==1:
            for key,branch in graph[item.k].adjk.items():#calcula as derivadas relativa as barras de e a do próprio angulo
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dQfdt(graph,0,item.k)
                if  branch.para in var_t.keys():
                    H[i][var_t[branch.para]]=branch.dQfdt(graph,0,branch.para)
            for key,branch in graph[item.k].adjm.items():
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dQfdt(graph,1,item.k) 
                if  branch.de in var_t.keys():
                    H[i][var_t[branch.de]]=branch.dQfdt(graph,1,branch.de)
            if  graph[item.k].bar.type!=0:
                H[i][var_t[item.k]]=soma1
            soma1=0
            for key,branch in graph[item.k].adjk.items(): ##derivadas modulo de tensão
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dQfdV(graph,0,item.k)
                if  branch.para in var_v.keys():
                    H[i][var_v[branch.para]+n_teta]=branch.dQfdV(graph,0,branch.para)
            for key,branch in graph[item.k].adjm.items():
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dQfdV(graph,1,item.k) 
                if  branch.de in var_v.keys():
                    H[i][var_v[branch.de]+n_teta]=branch.dQfdV(graph,1,branch.de)
            if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                if graph[item.k].FlagBS==1:
                    soma2=soma2-2*graph[item.k].Bs*graph[item.k].V 
                H[i][var_v[item.k]+n_teta]=soma2
            soma2=0
        elif item.type==2:
            print("teste")
        
        elif item.type==3:
            print("teste")

        elif item.type==4:
            print("teste")
        i=i+1
    return H

def calc_dz(vecZ,graph):
    dz=np.zeros(len(vecZ))
    i=0
    for z in vecZ:
        dz[i]=z.dz(graph)
        i=i+1
    return dz
def new_X(graph,var_t,var_v,dx):
    n_teta=len(var_t)
    for key,item in var_t.items():
        graph[key].teta=graph[key].teta+dx[item]
    for key,item in var_v.items():
        graph[key].V=graph[key].V+dx[item+n_teta]


def load_flow(graph,prt=0,tol=1e-6):
    Vinici_lf(graph)
    [z,var_t,var_v]=create_z_x_loadflow(graph)
    dx=np.ones(len(var_t))
    it=0
    conv=0
    while(it <20):
        dz=calc_dz(z,graph)
        H=calc_H(z,var_t,var_v,graph)
        if(it==0 and prt):
            np.savetxt("H.csv",H,delimiter=",")
        A=sparse.csc_matrix(H, dtype=float)
        dx=sliang.spsolve(A,dz)
        new_X(graph,var_t,var_v,dx)
        it=it+1
        if (np.amax(np.abs(dx))<tol):
            conv=1
            txt="Convergiu em {:d} iteracoes".format(it)
            print(txt)
            prt_state(graph)
            break
    if(it==0 and prt):
        np.savetxt("Hfinal.csv",H,delimiter=",")
    return conv


def create_z_x(graph,dfDMED,ind_i):
    z=[]
    var_t={}
    var_v={}
    i=0
    j=0
    for item in graph:
        if item.bar.type==1 or item.bar.type==2:
            var_t[item.id]=i
            i=i+1
        var_v[item.id]=j
        j=j+1

    for idx,row in dfDMED.iterrows():
        if int(row["type"])==0 or int(row["type"])==1 or  int(row["type"])==4:
            mes=meas(ind_i[int(row["de"])],-1,int(row["type"]),row["zmed"],row["prec"])
        else:  
            mes=meas(ind_i[int(row["de"])],ind_i[int(row["para"])],int(row["type"]),row["zmed"],row["prec"])
        z.append(mes)

    return z,var_t,var_v
