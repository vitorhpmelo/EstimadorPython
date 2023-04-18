from classes import *
from SS import *
import numpy as np
import pandas as pd
from readfiles import *
import scipy.sparse.linalg as sliang 
import scipy.sparse as sparse 
import csv




def SS_WLS_linear(graph,dfDMED,ind_i):

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    ref=list(set(list(range(len(graph))))-set(var_t.keys()))[0]
    zcc_conv=[]
    zcc_facts=[]
    for item in z:
        if item.type==0 or item.type==2:
            if item.type==2:
                k=item.k
                m=item.m
                if (str(k)+"-"+str(m)) in graph[k].bFACTS_adjk.keys() or (str(m)+"-"+str(k)) in graph[k].bFACTS_adjm.keys():
                    zcc_facts.append(item)
                else:
                    zcc_conv.append(item)
            else:
                zcc_conv.append(item)

    Hcc=np.zeros((len(zcc_conv)+len(zcc_facts),len(graph)))

    i=0
    d_injzcc={}
    for item in zcc_conv:
        if item.type==0:
            soma=0
            k=item.k
            d_injzcc[k]=i
            for key,ram in graph[k].adjk.items():
                if ram.type!=3:
                    Hcc[i][ram.para]=-1/(ram.x)
                    soma=soma+1/(ram.x)
            for key,ram in graph[k].adjm.items():
                if ram.type!=3:
                    Hcc[i][ram.de]=-1/(ram.x)
                    soma=soma+1/(ram.x)
            Hcc[i][k]=soma
        elif item.type==2:
            k=item.k
            m=item.m
            if str(k)+"-"+str(m) in graph[k].adjk.keys():
                ram=graph[k].adjk[str(k)+"-"+str(m)]
                Hcc[i][k]=1/(ram.x)
                Hcc[i][m]=-1/(ram.x)
            elif str(m)+"-"+str(k) in graph[k].adjm.keys():
                ram=graph[k].adjm[str(m)+"-"+str(k)]
                Hcc[i][k]=1/(ram.x)
                Hcc[i][m]=-1/(ram.x)
        i=i+1


    Hccx=np.zeros((len(zcc_conv)+len(zcc_facts),len(var_x)))

    for item in zcc_facts:
        k=item.k
        m=item.m
        Hcc[i][k]=1/item.val
        Hcc[i][m]=-1/item.val
        if str(k)+"-"+str(m) in var_x.keys():
            Hccx[i][var_x[str(k)+"-"+str(m)]]=-1
        elif str(m)+"-"+str(k) in var_x.keys():
            Hccx[i][var_x[str(m)+"-"+str(k)]]=-1
        i=i+1


    b=np.zeros(len(zcc_conv)+len(zcc_facts))
    i=0
    for item in zcc_conv:
        b[i]=item.val
        i=i+1        
    for item in zcc_facts:
        k=item.k
        m=item.m    
        if k in d_injzcc.keys():
            b[d_injzcc[k]]=b[d_injzcc[k]]-item.val
        if m in d_injzcc.keys():
            b[d_injzcc[m]]=b[d_injzcc[m]]+item.val
    
    W=create_W_cc(b,zcc_conv+zcc_facts,flag_ones=0)
    Hcc=np.concatenate((Hcc,Hccx),axis=1)

    H=np.delete(Hcc,ref,1)
    Gcc=np.matmul(H.T,np.matmul(W,H))
    b=np.matmul(np.matmul(H.T,W),b)
    x=np.linalg.solve(Gcc,b)
    return x,Hcc,var_t,var_v,var_x


def Vinici(graph,flatStart=0,dfDMED=[],ind_i=[]):
    '''
    Function to initate the voltages (state variables)
    If flagStart != 0 and != 1  with flat start (i.e. all the voltage modules equal to one and angles equal to 0)
    If flagStart == 0 the voltages from the DBAR
    If flagStart == 1 the voltages from the DBAR only in the ref bus
    @param: graph list of instances of the node class with all the information about the network
    @param: flagStart: 0 if the voltages should be initated with the values from the DBAR and different from 0 if they should initate with flat start  

    '''
    idxref=0
    for no in graph:
        if no.bar.type==0:
            idxref=no.id
            break

    if flatStart==0:
        for no in graph:
            no.V=no.bar.V
            no.teta=no.bar.teta
    elif flatStart==1:
        for no in graph:
            no.V=1
            no.teta=graph[idxref].bar.teta
    elif flatStart==2:
        for no in graph:
            no.V=no.bar.V
            no.teta=0
    elif flatStart==3:
        for no in graph:
            no.V=1+np.random.uniform(low=0,high=0.1)
            no.teta=0
    elif flatStart==4:
        for no in graph:
            if len(no.bFACTS_adjk)>0:
                no.V=1.1    
            else:
                no.V=1.0   
            no.teta=0
    elif flatStart==5:
        [x,H,var_t,var_v,var_x]=SS_WLS_linear(graph,dfDMED,ind_i)
        for no in graph:
            k=no.id
            if no.bar.type!=0:
                i=var_t[k]
                no.teta=x[i]
            no.V=1
        for key,i in var_x.items():
            k=int(key.split("-")[0])
            if key in graph[k].adjk.keys():
                graph[k].adjk[key].xtcsc=x[i+len(var_t)]
                graph[k].adjk[key].AttY()
    else:
        for no in graph:
            no.V=1
            no.teta=0

def Vinici_lf(graph,useDBAR=1,var_x=[],var_t=[],z=[]):
    '''
    Function to initate the voltages (state variables) for the load flow, 
    PQ buses recive 1 for the voltage module and 0 for the angle,
    PV recive the V from the DBAR for the module
    slack initate with the voltage from the DB 
    @param: graph list of instances of the node class with all the information about the network
    '''
    tetaini=0
    for no in graph:
        if no.bar.type == 0:
            tetaini=no.bar.teta
            break

    if useDBAR!=-1:
        for no in graph:
            if no.bar.type == 1 or no.bar.type == 0:
                if useDBAR==0:
                    no.V=1
                    no.teta=0
                elif useDBAR==1:
                    no.V=no.bar.V
                    if no.bar.type == 0:
                        no.teta=no.bar.teta
                    else:
                        no.teta=0+tetaini  
            else:
                no.V=1
                no.teta=0+tetaini 
    else:
        [x,H]=load_flow_FACTS_cc(z,graph,var_x,var_t)
        for no in graph:
            k=no.id
            if no.bar.type!=0:
                i=var_t[k]
                no.teta=x[i]+tetaini
            else:
                no.teta=tetaini
            if no.bar.type == 1 or no.bar.type == 0:
                no.V=no.bar.V
            else:
                no.V=1
        for key,i in var_x.items():
            k=int(key.split("-")[0])
            if key in graph[k].adjk.keys():
                graph[k].adjk[key].xtcsc=x[i+len(var_t)]
                graph[k].adjk[key].AttY()




def Vinici_DBAR(graph):
    '''
    Function to initate the voltages (state variables) for the load flow, 
    PQ buses recive 1 for the voltage module and 0 for the angle,
    PV recive the V from the DBAR for the module
    slack initate with the voltage from the DB 
    @param: graph list of instances of the node class with all the information about the network
    '''
    for no in graph:
        no.V=no.bar.V
        no.teta=round(no.bar.teta,1)

def FACTSini(graph,useDFACTS=1):
    """
    Function to initialize FACTS devices
    """
    if useDFACTS==0:
        for no in graph:
            if no.FlagTCSC==1:
                for  key in no.bFACTS_adjk.keys():
                    no.bFACTS_adjk[key].xtcsc=0.0001
                    no.bFACTS_adjk[key].AttY()
                    
    if useDFACTS==1:
        for no in graph:
            if no.FlagTCSC==1:
                for  key in no.bFACTS_adjk.keys():
                    no.bFACTS_adjk[key].xtcsc=no.bFACTS_adjk[key].xtcsc_ini
                    no.bFACTS_adjk[key].AttY()



def PowerFlows(ram,graph,print=0):
    """
    Funtion to calculate the power flows acros the all the network branches and return them into the dpf (for active) and dqf (for reactive)
    dictionary 
    """
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

def create_z_x_loadflow_TCSC(graph):
    zPf=[]
    var_x={}
    i=0
    for no in graph:
        if no.FlagTCSC==1 and len(no.bFACTS_adjk.keys())>0:
            for key,item in no.bFACTS_adjk.items():
                mes=meas(item.de,item.para,2,item.Pfesp,1)
                zPf.append(mes)
                var_x[str(item.de)+"-"+str(item.para)]=i
                i=i+1
    return zPf,var_x

def create_x_TCSC(graph):
    var_x={}
    i=0
    for no in graph:
        if no.FlagTCSC==1 and len(no.bFACTS_adjk.keys())>0:
            for key,item in no.bFACTS_adjk.items():
                var_x[str(item.de)+"-"+str(item.para)]=i
                i=i+1
    return var_x

def calc_H_fp(z,var_t,var_v,graph,H):
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
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].adjk[km].dPfdt(graph,0,k)
                if k in var_v.keys():
                    H[i][var_v[k]+n_teta]= graph[k].adjk[km].dPfdV(graph,0,k)
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].adjk[km].dPfdt(graph,0,m)
                if m in var_v.keys():
                    H[i][var_v[m]+n_teta]= graph[k].adjk[km].dPfdV(graph,0,m)
            elif mk in graph[k].adjm.keys():
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].adjm[mk].dPfdt(graph,1,k)
                H[i][var_v[k]+n_teta]= graph[k].adjm[mk].dPfdV(graph,1,k)
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].adjm[mk].dPfdt(graph,1,m)
                H[i][var_v[m]+n_teta]= graph[k].adjm[mk].dPfdV(graph,1,m)
            else:
                print("erro ao calcular fluxo na Jacobiana, medida Fluxo de P {:d}-{:d}".format(graph[k].id,graph[m].id))
                exit(1)
        elif item.type==3:
                k=item.k
                m=item.m
                km=str(k)+"-"+str(m)
                mk=str(m)+"-"+str(k)
                if km in graph[k].adjk.keys():
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].adjk[km].dQfdt(graph,0,k)
                    H[i][var_v[k]+n_teta]= graph[k].adjk[km].dQfdV(graph,0,k)
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].adjk[km].dQfdt(graph,0,m)
                    H[i][var_v[m]+n_teta]= graph[k].adjk[km].dQfdV(graph,0,m)
                elif mk in graph[k].adjm.keys():
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].adjm[mk].dQfdt(graph,1,k)
                    H[i][var_v[k]+n_teta]= graph[k].adjm[mk].dQfdV(graph,1,k)
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].adjm[mk].dQfdt(graph,1,m)
                    H[i][var_v[m]+n_teta]= graph[k].adjm[mk].dQfdV(graph,1,m)
                else:
                    print("erro ao calcular fluxo na Jacobiana, medida Fluxo de P {:d}-{:d}".format(graph[k].id,graph[m].id))
                    exit(1)
        elif item.type==4:
                k=item.k
                H[i][var_v[k]+n_teta]=1
        i=i+1


def calc_H_fp_TCSC(z,var_x,graph,H):
    i=0

    for item in z:
        if item.type==0:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dPfdx(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dPfdx(graph,1)     
        elif item.type==1:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dQfdx(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dQfdx(graph,1)    
        elif item.type==2:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                H[i][var_x[km]]= graph[k].adjk[km].dPfdx(graph,0)
            elif mk in graph[k].adjm.keys():
                H[i][var_x[mk]]= graph[k].adjm[mk].dPfdx(graph,1)
        elif item.type==3:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                H[i][var_x[km]]= graph[k].adjk[km].dQfdx(graph,0)
            elif mk in graph[k].adjm.keys():
                H[i][var_x[mk]]= graph[k].adjm[mk].dQfdx(graph,1)
            else:
                print("erro ao calcular fluxo na Jacobiana, medida Fluxo deQ {:d}-{:d}".format(graph[k].id,graph[m].id))
                exit(1)
        elif item.type==4:
                for key in var_x.keys():
                    H[i][var_x[key]]=0 
        i=i+1


def calc_H_fp_TCSC_B(z,var_x,graph,H):
    i=0

    for item in z:
        if item.type==0:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dPfdB(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dPfdB(graph,1)     
        elif item.type==1:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dQfdB(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dQfdB(graph,1)    
        elif item.type==2:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                H[i][var_x[km]]= graph[k].adjk[km].dPfdB(graph,0)
            elif mk in graph[k].adjm.keys():
                H[i][var_x[mk]]= graph[k].adjm[mk].dPfdB(graph,1)
        elif item.type==3:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                H[i][var_x[km]]= graph[k].adjk[km].dQfdB(graph,0)
            elif mk in graph[k].adjm.keys():
                H[i][var_x[mk]]= graph[k].adjm[mk].dQfdB(graph,1)
            else:
                print("erro ao calcular fluxo na Jacobiana, medida Fluxo deQ {:d}-{:d}".format(graph[k].id,graph[m].id))
                exit(1)
        elif item.type==4:
                for key in var_x.keys():
                    H[i][var_x[key]]=0 
        i=i+1


def calc_H_EE_TCSC(z,var_x,graph,H):
    i=0

    for item in z:
        if item.type==0:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dPfdx(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dPfdx(graph,1)     
        elif item.type==1:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dQfdx(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dQfdx(graph,1)    
        elif item.type==2:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                if  graph[k].adjk[km].type==3:
                    H[i][var_x[km]]= graph[k].adjk[km].dPfdx(graph,0)
            elif mk in graph[k].adjm.keys():
                if graph[k].adjm[mk].type==3:
                    H[i][var_x[mk]]= graph[k].adjm[mk].dPfdx(graph,1)
        elif item.type==3:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                if  graph[k].adjk[km].type==3:
                    H[i][var_x[km]]= graph[k].adjk[km].dQfdx(graph,0)
            elif mk in graph[k].adjm.keys():
                if graph[k].adjm[mk].type==3:
                    H[i][var_x[mk]]= graph[k].adjm[mk].dQfdx(graph,1)    
            else:
                print("erro ao calcular fluxo na Jacobiana, medida Fluxo deQ {:d}-{:d}".format(graph[k].id,graph[m].id))
                exit(1)
        elif item.type==4:
                for key in var_x.keys():
                    H[i][var_x[key]]=0 
        i=i+1



def calc_H_EE_TCSC_B(z,var_x,graph,H):
    i=0

    for item in z:
        if item.type==0:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dPfdB(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dPfdB(graph,1)     
        elif item.type==1:
            k=item.k
            for key in set(graph[k].adjk.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjk[key].dQfdB(graph,0)
            for key in set(graph[k].adjm.keys()).intersection(set(var_x.keys())):
                H[i][var_x[key]]=graph[k].adjm[key].dQfdB(graph,1)    
        elif item.type==2:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                if  graph[k].adjk[km].type==3:
                    H[i][var_x[km]]= graph[k].adjk[km].dPfdB(graph,0)
            elif mk in graph[k].adjm.keys():
                if graph[k].adjm[mk].type==3:
                    H[i][var_x[mk]]= graph[k].adjm[mk].dPfdB(graph,1)
        elif item.type==3:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                if  graph[k].adjk[km].type==3:
                    H[i][var_x[km]]= graph[k].adjk[km].dQfdB(graph,0)
            elif mk in graph[k].adjm.keys():
                if graph[k].adjm[mk].type==3:
                    H[i][var_x[mk]]= graph[k].adjm[mk].dQfdB(graph,1)    
            else:
                print("erro ao calcular fluxo na Jacobiana, medida Fluxo deQ {:d}-{:d}".format(graph[k].id,graph[m].id))
                exit(1)
        elif item.type==4:
                for key in var_x.keys():
                    H[i][var_x[key]]=0 
        i=i+1

def calc_H_EE(z,var_t,var_v,graph,H):
    #refazer
    i=0
    n_teta=len(var_t)
    bar_v=var_v.keys()
    bar_t=var_t.keys()
    for item in z:
        soma1=0
        soma2=0
        if item.type==0:
            for key,branch in graph[item.k].adjk.items():#calcula as derivadas relativa as barras de e a do próprio angulo
                if  item.k in bar_t: #checa se o angulo daquela barra é variável
                    soma1=soma1+branch.dPfdt(graph,0,item.k) #soma a derivada de cada fluxo incidente a barra de
                if  branch.para in bar_t:
                    H[i][var_t[branch.para]]=branch.dPfdt(graph,0,branch.para)  # calcula a derivada daquela barra em relação a seguinte
            for key,branch in graph[item.k].adjm.items(): #calcula as derivadas em relação aos fluxo em relação aos nos que aquela barra é para
                if  item.k in bar_t: #checa se o angulo daquela barra é variável
                    soma1=soma1+branch.dPfdt(graph,1,item.k) #soma a derivada de cada fluxo incidente a barra de
                if  branch.de in bar_t:#checa se o angulo daquela barra é variável
                    H[i][var_t[branch.de]]=branch.dPfdt(graph,1,branch.de) 
            if  item.k in bar_t:
                H[i][var_t[item.k]]=soma1
            soma1=0
            for key,branch in graph[item.k].adjk.items(): ##derivadas modulo de tensão
                if item.k in bar_v:
                    soma2=soma2+branch.dPfdV(graph,0,item.k)
                if  branch.de in var_v.keys():
                    H[i][var_v[branch.para]+n_teta]=branch.dPfdV(graph,0,branch.para)
            for key,branch in graph[item.k].adjm.items():
                if item.k in bar_v:
                    soma2=soma2+branch.dPfdV(graph,1,item.k) 
                if  branch.de in var_v.keys():
                    H[i][var_v[branch.de]+n_teta]=branch.dPfdV(graph,1,branch.de)
            if item.k in bar_v:
                H[i][var_v[item.k]+n_teta]=soma2
            soma2=0
        elif item.type==1:
            for key,branch in graph[item.k].adjk.items():#calcula as derivadas relativa as barras de e a do próprio angulo
                if  item.k in bar_t:
                    soma1=soma1+branch.dQfdt(graph,0,item.k)
                if  branch.para in bar_t:
                    H[i][var_t[branch.para]]=branch.dQfdt(graph,0,branch.para)
            for key,branch in graph[item.k].adjm.items():
                if  item.k in bar_t:
                    soma1=soma1+branch.dQfdt(graph,1,item.k) 
                if  branch.de in bar_t:
                    H[i][var_t[branch.de]]=branch.dQfdt(graph,1,branch.de)
            if  item.k in bar_t:
                H[i][var_t[item.k]]=soma1
            soma1=0
            for key,branch in graph[item.k].adjk.items(): ##derivadas modulo de tensão
                if item.k in bar_v:
                    soma2=soma2+branch.dQfdV(graph,0,item.k)
                if  branch.para in bar_v:
                    H[i][var_v[branch.para]+n_teta]=branch.dQfdV(graph,0,branch.para)
            for key,branch in graph[item.k].adjm.items():
                if item.k in bar_v:
                    soma2=soma2+branch.dQfdV(graph,1,item.k) 
                if  branch.de in bar_v:
                    H[i][var_v[branch.de]+n_teta]=branch.dQfdV(graph,1,branch.de)
            if graph[item.k].FlagBS==1:
                soma2=soma2-2*graph[item.k].Bs*graph[item.k].V 
            H[i][var_v[item.k]+n_teta]=soma2
            soma2=0
        elif item.type==2:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].adjk.keys():
                if k in bar_t:
                    H[i][var_t[k]]= graph[k].adjk[km].dPfdt(graph,0,k)
                H[i][var_v[k]+n_teta]= graph[k].adjk[km].dPfdV(graph,0,k)
                if m in bar_t:
                    H[i][var_t[m]]= graph[k].adjk[km].dPfdt(graph,0,m)
                H[i][var_v[m]+n_teta]= graph[k].adjk[km].dPfdV(graph,0,m)
            elif mk in graph[k].adjm.keys():
                if k in bar_t:
                    H[i][var_t[k]]= graph[k].adjm[mk].dPfdt(graph,1,k)
                H[i][var_v[k]+n_teta]= graph[k].adjm[mk].dPfdV(graph,1,k)
                if m in bar_t:
                    H[i][var_t[m]]= graph[k].adjm[mk].dPfdt(graph,1,m)
                H[i][var_v[m]+n_teta]= graph[k].adjm[mk].dPfdV(graph,1,m)
            else:
                print("erro ao calcular fluxo na Jacobiana, medida Fluxo de P {:d}-{:d}".format(graph[k].id,graph[m].id))
                exit(1)
        elif item.type==3:
                k=item.k
                m=item.m
                km=str(k)+"-"+str(m)
                mk=str(m)+"-"+str(k)
                if km in graph[k].adjk.keys():
                    if k in bar_t:
                        H[i][var_t[k]]= graph[k].adjk[km].dQfdt(graph,0,k)
                    H[i][var_v[k]+n_teta]= graph[k].adjk[km].dQfdV(graph,0,k)
                    if m in bar_t:
                        H[i][var_t[m]]= graph[k].adjk[km].dQfdt(graph,0,m)
                    H[i][var_v[m]+n_teta]= graph[k].adjk[km].dQfdV(graph,0,m)
                elif mk in graph[k].adjm.keys():
                    if k in bar_t:
                        H[i][var_t[k]]= graph[k].adjm[mk].dQfdt(graph,1,k)
                    H[i][var_v[k]+n_teta]= graph[k].adjm[mk].dQfdV(graph,1,k)
                    if m in bar_t:
                        H[i][var_t[m]]= graph[k].adjm[mk].dQfdt(graph,1,m)
                    H[i][var_v[m]+n_teta]= graph[k].adjm[mk].dQfdV(graph,1,m)
                else:
                    print("erro ao calcular fluxo na Jacobiana, medida Fluxo de P {:d}-{:d}".format(graph[k].id,graph[m].id))
                    exit(1)
        elif item.type==4:
                k=item.k
                H[i][var_v[k]+n_teta]=1
        i=i+1


def calc_dz(vecZ,graph,dz):
    i=0
    for z in vecZ:
        dz[i]=z.dz(graph)
        i=i+1


def calc_cx(vecc,graph,cx):
    i=0
    for c in vecc:
        cx[i]=c.cx(graph)
        i=i+1


def new_X(graph,var_t,var_v,dx):
    n_teta=len(var_t)
    for key,item in var_t.items():
        graph[key].teta=graph[key].teta+dx[item]
    for key,item in var_v.items():
        graph[key].V=graph[key].V+dx[item+n_teta]

def new_X_TCSCC(graph,nvars,var_x,dx):
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        graph[k].adjk[key].xtcsc=graph[k].adjk[key].xtcsc+dx[item+nvars]
        graph[k].adjk[key].AttY()


def new_X_TCSCC_B(graph,nvars,var_x,dx):
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        graph[k].adjk[key].xtcsc=1/(1/(graph[k].adjk[key].xtcsc)+dx[item+nvars])
        graph[k].adjk[key].AttY()





def load_flow_FACTS(graph,prt=0,tol=1e-6,inici=-1,itmax=20):
    """
    Function to run load flow with FACTS devices (only TCSC implemented yet)
    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance of the load flow calculation
    @param inici intialization method of the variables if -1 it uses the DC power flow to intialize the angles and the X of the TCSC, 
    if it is 1 other value it uses DBAR for the PV voltage magnitudes and if it is 0, it initalizes with flat start
    """

    
    zPf,var_x = create_z_x_loadflow_TCSC(graph)#create z and dinctionary with the variables for the FACTS devices
    [z,var_t,var_v]=create_z_x_loadflow(graph)#create z and var_v and var_t for the traditional load flow
    z=z+zPf
    FACTSini(graph,useDFACTS=1)
    Vinici_lf(graph,useDBAR=inici,var_t=var_t,var_x=var_x,z=z)
    dz=np.zeros(len(z))
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    Hx=np.zeros((len(z),len(var_x)))
    it=0
    conv=0
    lstdx=[]
    lstdz=[]

    while it<itmax:
        calc_dz(z,graph,dz)
        calc_H_fp(z,var_t,var_v,graph,H)
        calc_H_fp_TCSC(z,var_x,graph,Hx)
        HTCSC=np.concatenate((H,Hx),axis=1)
        A=sparse.csc_matrix(HTCSC, dtype=float)
        dx=sliang.spsolve(A,dz)
        new_X(graph,var_t,var_v,dx)
        new_X_TCSCC(graph,len(var_t)+len(var_v),var_x,dx)
        maxdx=np.max(np.abs(dx))
        maxdz=np.max(np.abs(dz))
        lstdx.append(maxdx)
        lstdz.append(maxdz)
        if maxdx< tol and maxdz < tol:
            print("convergiu em {} itereacoes".format(it))
            prt_state(graph)
            conv=1
            break
        it=it+1

    if prt==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv.csv', index=False)
    return conv


def load_flow_FACTS_2(graph,prt=0,tol=1e-6,inici=-1,itmax=20):
    '''
    Function to run load flow with FACTS devices (only TCSC implemented yet)
    run the load flow but it uses uses the B as the state variable
    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance of the load flow calculation
    @param inici intialization method of the variables if -1 it uses the DC power flow to intialize the angles and the X of the TCSC, 
    if it is 1 other value it uses DBAR for the PV voltage magnitudes and if it is 0, it initalizes with flat start
    
    '''

    zPf,var_x = create_z_x_loadflow_TCSC(graph)
    [z,var_t,var_v]=create_z_x_loadflow(graph)
    z=z+zPf
    FACTSini(graph,useDFACTS=1)
    Vinici_lf(graph,useDBAR=inici,var_t=var_t,var_x=var_x,z=z)
    dz=np.zeros(len(z))
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    Hx=np.zeros((len(z),len(var_x)))
    it=0
    conv=0
    lstdx=[]
    lstdz=[]
    while it<20:
        calc_dz(z,graph,dz)
        calc_H_fp(z,var_t,var_v,graph,H)
        calc_H_fp_TCSC_B(z,var_x,graph,Hx)
        HTCSC=np.concatenate((H,Hx),axis=1)
        A=sparse.csc_matrix(HTCSC, dtype=float)
        dx=sliang.spsolve(A,dz)
        new_X(graph,var_t,var_v,dx)
        new_X_TCSCC_B(graph,len(var_t)+len(var_v),var_x,dx)
        maxdx=np.max(np.abs(dx))
        maxdz=np.max(np.abs(dz))
        lstdx.append(maxdx)
        lstdz.append(maxdz)
        if np.max(np.abs(dx))< tol and np.max(np.abs(dz)) < tol:
            print("convergiu em {} itereacoes".format(it))
            prt_state(graph)
            conv=1
            break
        it=it+1

    if prt==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)
        # Save the DataFrame to a CSV file
        df.to_csv('conv.csv', index=False)

    return conv


def load_flow(graph,prt=0,tol=1e-6):
    Vinici_lf(graph)
    
    [z,var_t,var_v]=create_z_x_loadflow(graph)
    
    dx=np.ones(len(var_t))
    dz=np.zeros(len(z))
    it=0
    conv=0
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    while(it <20):
        calc_dz(z,graph,dz)
        calc_H_fp(z,var_t,var_v,graph,H)
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



def create_z_c_x_LGI(graph,dfDMED,ind_i):
    z=[]
    c=[]
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
        if (int(row["type"])==0 or int(row["type"])==1) and row["zmed"]==0:
            c.append(mes)
        else:
            z.append(mes)

    return z,c,var_t,var_v

def create_W(z,prec_virtual=1e-4,flag_ones=0):

    if flag_ones==0:
        W=np.zeros((len(z),len(z)))
        i=0
        for item in z:
            if np.abs(item.val)<1e-6:
                W[i][i]=1/(prec_virtual**2)
            else:
                W[i][i]=1/(item.sigma**2)
            i=i+1
    else:
        W=np.eye(len(z))
    return W


def create_W_cc(b,z,prec_virtual=1e-4,flag_ones=0):

    if flag_ones==0:
        W=np.zeros((len(z),len(z)))
        i=0
        for item in b:
            if np.abs(item)<1e-6:
                W[i][i]=1/(prec_virtual**2)
            else:
                sigma=np.abs(item)*z[i].prec/3
                W[i][i]=1/(sigma**2)
            i=i+1
    else:
        W=np.eye(len(z))
    return W


def calcYbus(graph,ram):
    YBus=np.zeros([len(graph),len(graph)],dtype=complex)

    for key,item in ram.items():
        i=item.de
        j=item.para
        YBus[i][j]=YBus[i][j]+item.Y[0][1]
        YBus[j][i]=YBus[j][i]+item.Y[1][0]
        YBus[i][i]=YBus[i][i]+item.Y[0][0]
        YBus[j][j]=YBus[j][j]+item.Y[1][1]

    for no in graph:
        if no.FlagBS==1:
            i=no.id
            YBus[i][i]=YBus[i][i]+complex(0,no.Bs)

    np.savetxt("Ybus.csv",YBus)


def load_flow_FACTS_cc(z,graph,var_x,var_t):
    ref=list(set(list(range(len(graph))))-set(var_t.keys()))[0]
    zcc_conv=[]
    zcc_facts=[]
    for item in z:
        if item.type==0 or item.type==2:
            if item.type==2:
                k=item.k
                m=item.m
                if (str(k)+"-"+str(m)) in graph[k].bFACTS_adjk.keys() or (str(m)+"-"+str(k)) in graph[k].bFACTS_adjm.keys():
                    zcc_facts.append(item)
                else:
                    zcc_conv.append(item)
            else:
                zcc_conv.append(item)

    Hcc=np.zeros((len(zcc_conv)+len(zcc_facts),len(graph)))

    i=0
    d_injzcc={}
    for item in zcc_conv:
        if item.type==0:
            soma=0
            k=item.k
            d_injzcc[k]=i
            for key,ram in graph[k].adjk.items():
                if ram.type!=3:
                    Hcc[i][ram.para]=-1/(ram.x)
                    soma=soma+1/(ram.x)
            for key,ram in graph[k].adjm.items():
                if ram.type!=3:
                    Hcc[i][ram.de]=-1/(ram.x)
                    soma=soma+1/(ram.x)
            Hcc[i][k]=soma
        elif item.type==2:
            k=item.k
            m=item.m
            Hcc[i][k]=1/(ram.x)
            Hcc[i][m]=-1/(ram.x)
        i=i+1

    Hccx=np.zeros((len(zcc_conv)+len(zcc_facts),len(var_x)))
    for item in zcc_facts:
        k=item.k
        m=item.m
        Hcc[i][k]=1/item.val
        Hcc[i][m]=-1/item.val
        if str(k)+"-"+str(m) in var_x.keys():
            Hccx[i][var_x[str(k)+"-"+str(m)]]=-1
        elif str(m)+"-"+str(k) in var_x.keys():
            Hccx[i][var_x[str(k)+"-"+str(m)]]=-1
        i=i+1

    b=np.zeros(len(zcc_conv)+len(zcc_facts))
    i=0
    for item in zcc_conv:
        b[i]=item.val
        i=i+1
    for item in zcc_facts:
        k=item.k
        m=item.m    
        if k in d_injzcc.keys():
            b[d_injzcc[k]]=b[d_injzcc[k]]-item.val
        if m in d_injzcc.keys():
            b[d_injzcc[m]]=b[d_injzcc[m]]+item.val
    A=np.concatenate((Hcc,Hccx),axis=1)

    A_red=np.delete(A,ref,1)

    x=np.linalg.solve(A_red,b)
    return x,A