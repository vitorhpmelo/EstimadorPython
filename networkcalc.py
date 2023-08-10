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
    for no in graph:
        if no.bar.type == 0:
            tetaini=no.bar.teta
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
            no.teta=tetaini
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
    elif flatStart==5:
        [x,H,var_t,var_v,var_x]=SS_WLS_linear(graph,dfDMED,ind_i)
        for no in graph:
            k=no.id
            if no.bar.type!=0:
                i=var_t[k]
                no.teta=x[i]+tetaini
            elif no.bar.type==0:
                no.teta=tetaini
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
                if  no.FlagSVC==True:
                    no.V=no.bar.V 
                elif (no.FlagSVC==False) and (no.FlagUPFC==True) and (len(no.bUFPC_adjk.items())>0):
                    no.V=no.bar.V 
                else:
                    no.V=1
                no.teta=0+tetaini 
                    
        for key in var_x.keys():
            key=key.split("-")
            k=int(key[1])
            graph[k].V=graph[k].V+0.1

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
                if  no.FlagSVC==True:
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
        no.teta=no.bar.teta

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
            if no.FlagSVC==True:
                no.SVC.BSVC=0.10
                no.SVC.attYk()    
            if no.FlagUPFC==1:
                for  key in no.bUFPC_adjk.keys():
                    no.bUFPC_adjk[key].Vse=0.02
                    no.bUFPC_adjk[key].Vsh=1.00
                    no.bUFPC_adjk[key].t_se=-90/np.pi()
                    no.bUFPC_adjk[key].t_sh=0
    if useDFACTS==1:
        for no in graph:
            if no.FlagTCSC==1:
                for  key in no.bFACTS_adjk.keys():
                    no.bFACTS_adjk[key].xtcsc=no.bFACTS_adjk[key].xtcsc_ini
                    no.bFACTS_adjk[key].AttY()
            if no.FlagSVC==True:
                no.SVC.BSVC=no.SVC.Bini
                no.SVC.attYk() 
            if no.FlagUPFC==True:
                for key, item in no.bUFPC_adjk.items():
                    if item.mode==1:
                        no.V=item.Vp 
            if no.FlagUPFC==1:
                for  key in no.bUFPC_adjk.keys():
                    no.bUFPC_adjk[key].Vse=no.bUFPC_adjk[key].Vse_ini
                    no.bUFPC_adjk[key].Vsh=no.bUFPC_adjk[key].Vsh_ini
                    no.bUFPC_adjk[key].t_se=no.bUFPC_adjk[key].t_se_ini
                    no.bUFPC_adjk[key].t_sh=no.bUFPC_adjk[key].t_sh_ini



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


def create_x_loadflow_SVC(graph,var_v):
    """
    Auxiliar funcition for the load flow routine, creates the dictionary with the variables from the SVC and REMOVES the voltage magnitudes of the buses with SVC.

    """

    
    var_svc={}
    i=0
    for no in graph:
        if no.FlagSVC==True:
            var_svc[no.id]=i
            i=i+1
    
    i=0
    keys=list(var_v.keys())
    for key in keys:
        if key in var_svc.keys():
            del var_v[key]
        else:
            var_v[key]=i
            i=i+1

    return var_svc


def create_x_loadflow_UPFC(graph,var_v):
    """
    Auxiliar funcition for the load flow routine, creates the dictionary with the variables from the SVC and REMOVES the voltage magnitudes of the buses with SVC.

    """

    
    var_UPFC={}
    var_UPFC_vsh={}
    i=0
    i_vsh=0
    zPf=[]
    zQf=[]
    p_controlled=[]
    for no in graph:
        if no.FlagUPFC==1 and len(no.bUFPC_adjk.keys())>0:
            for key,item in no.bUFPC_adjk.items():
                mes=meas(item.s,item.p,2,item.Psp_set,1)
                zPf.append(mes)
                mes=meas(item.s,item.p,3,item.Qsp_set,1)
                zQf.append(mes)
                var_UPFC[str(item.p)+"-"+str(item.s)]=i
                i=i+1
                if item.mode==1:
                    var_UPFC_vsh[str(item.p)+"-"+str(item.s)]=i_vsh
                    i_vsh=i_vsh+1
                    p_controlled.append(item.p)


    i=0
    keys=list(var_v.keys())
    for key in keys:
        if key in p_controlled:
            del var_v[key]
        else:
            var_v[key]=i
            i=i+1

    return [zPf+zQf,var_UPFC,var_UPFC_vsh]

def create_x_SVC(graph):
    """
    Creates the dictionary for the SVC variables
    """
    
    var_svc={}
    i=0
    for no in graph:
        if no.FlagSVC==True:
            var_svc[no.id]=i
            i=i+1
    
    return var_svc


def create_x_TCSC(graph):
    var_x={}
    i=0
    for no in graph:
        if no.FlagTCSC==1 and len(no.bFACTS_adjk.keys())>0:
            for key,item in no.bFACTS_adjk.items():
                var_x[str(item.de)+"-"+str(item.para)]=i
                i=i+1
    return var_x


def create_c_x_UPFC(graph):
    """
    Auxiliar funcition for the load flow routine, creates the dictionary with the variables for the UPFC 
    for the State Estimator
    """
    var_UPFC={}
    i=0
    c_upfc=[]
    for no in graph:
        if no.FlagUPFC==1 and len(no.bUFPC_adjk.keys())>0:
            for key,item in no.bUFPC_adjk.items():
                var_UPFC[str(item.p)+"-"+str(item.s)]=i
                i=i+1


    c_upfc=np.zeros(len(var_UPFC))


    return var_UPFC,c_upfc


def calc_H_fp(z,var_t,var_v,graph,H):
    i=0
    n_teta=len(var_t)
    for item in z:
        soma1=0
        soma2=0
        if item.type==0:
            #-------------------ramos de branchs normais e TCSC----------------------------------------#
            for key,branch in graph[item.k].adjk.items():# o branch entra com k-m e barra k é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dPfdt(graph,0,item.k) # cacula dPkm/dtk
                if  branch.para in var_t.keys():
                    H[i][var_t[branch.para]]=branch.dPfdt(graph,0,branch.para) #caclula dPkm/dtm para teta m na jacobiana
            for key,branch in graph[item.k].adjm.items(): # o branch entra com k-m e barra m é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dPfdt(graph,1,item.k)  # calcula dpmk/dm
                if  branch.de in var_t.keys():
                    H[i][var_t[branch.de]]=branch.dPfdt(graph,1,branch.de) #faz calcula dPmk/dk
            #-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+upfc.dPpsdtp(graph) # cacula dPps/dtp dPfdt(graph,0,de)
                if upfc.s in var_t.keys():
                    H[i][var_t[upfc.s]]=upfc.dPpsdts(graph) # cacula dPps/dts  dPfdt(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+upfc.dPspdts(graph) # calcula dPsp/dts dfPf(graph,1,para)
                if upfc.p in var_t.keys():
                    H[i][var_t[upfc.p]]=upfc.dPspdtp(graph) #cacula dPsp/dtp  dfPf(graph,1,de)
            if  graph[item.k].bar.type!=0:
                H[i][var_t[item.k]]=soma1
            soma1=0
            #-------------------ramos de branchs normais e TCSC----------------------------------------#
            for key,branch in graph[item.k].adjk.items(): # o branch entra com k-m e barra k é a variável
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dPfdV(graph,0,item.k) # calcula dPkm/dVk
                if  branch.para in var_v.keys():
                    H[i][var_v[branch.para]+n_teta]=branch.dPfdV(graph,0,branch.para) # calcula dPkm/dVm
            for key,branch in graph[item.k].adjm.items():  # o branch entra com k-m e barra m é a variável
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dPfdV(graph,1,item.k) # Calcula dPmk/dm 
                if  branch.de in var_v.keys():
                    H[i][var_v[branch.de]+n_teta]=branch.dPfdV(graph,1,branch.de) # Calcula dPmk/dk
            #-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if  graph[item.k].bar.type!=0:
                    soma2=soma2+upfc.dPpsdVp(graph) # cacula dPps/dVp dPfdV(graph,0,de)
                if upfc.s in var_v.keys():
                    H[i][var_v[upfc.s]+n_teta]=upfc.dPpsdVs(graph) # cacula dPps/dVs  dPfdV(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  graph[item.k].bar.type!=0:
                    soma2=soma2+upfc.dPspdVs(graph) # calcula dPsp/dVs dfPf(graph,1,para)
                if upfc.p in var_v.keys():
                    H[i][var_v[upfc.p]+n_teta]=upfc.dPspdVp(graph) #cacula dPsp/dVp  dfPfdV(graph,1,de)          
            if graph[item.k].FlagSVC==1:
                soma2=soma2+ 2*graph[item.k].SVC.Gk*graph[item.k].V
            if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1 and (item.k in var_v.keys()):
                H[i][var_v[item.k]+n_teta]=soma2## Colocar as derivadas do shunt
            soma2=0
        elif item.type==1:
            #-------------------ramos de branchs normais e TCSC----------------------------------------#
            for key,branch in graph[item.k].adjk.items(): # o branch entra com k-m e barra k é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dQfdt(graph,0,item.k) # caclula dQkm/dtk
                if  branch.para in var_t.keys(): 
                    H[i][var_t[branch.para]]=branch.dQfdt(graph,0,branch.para) #caclula da dQkm/dtm
            for key,branch in graph[item.k].adjm.items(): # o branch entra com k-m e barra m é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+branch.dQfdt(graph,1,item.k)  #caclula dQmk/dtm
                if  branch.de in var_t.keys():
                    H[i][var_t[branch.de]]=branch.dQfdt(graph,1,branch.de) #cacula dQmk/dk
            #-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+upfc.dQpsdtp(graph) # cacula dQps/dtp dQfdt(graph,0,de)
                if upfc.s in var_t.keys():
                    H[i][var_t[upfc.s]]=upfc.dQpsdts(graph) # cacula dQps/dts  dQfdt(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  graph[item.k].bar.type!=0:
                    soma1=soma1+upfc.dQspdts(graph) # calcula dQsp/dts dfQf(graph,1,para)
                if upfc.p in var_t.keys():
                    H[i][var_t[upfc.p]]=upfc.dQspdtp(graph) #cacula dQsp/dtp  dfPf(graph,1,de)  
            if  graph[item.k].bar.type!=0:
                H[i][var_t[item.k]]=soma1
            soma1=0
            #-------------------ramos de branchs normais e TCSC----------------------------------------#
            for key,branch in graph[item.k].adjk.items(): ## o branch entra com p-s e barra p é a variável (derivadas modulo de tensão)
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dQfdV(graph,0,item.k) #caclula dQkm/dVk
                if  branch.para in var_v.keys():
                    H[i][var_v[branch.para]+n_teta]=branch.dQfdV(graph,0,branch.para) #cacula dQkm/dVm
            for key,branch in graph[item.k].adjm.items(): ## o branch entra com p-s e barra s é a variável (derivadas modulo de tensão)
                if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1:
                    soma2=soma2+branch.dQfdV(graph,1,item.k) # caclula dQmk/dVm
                if  branch.de in var_v.keys():
                    H[i][var_v[branch.de]+n_teta]=branch.dQfdV(graph,1,branch.de) #cacula dQmk/dVk
            #-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if  graph[item.k].bar.type!=0:
                    soma2=soma2+upfc.dQpsdVp(graph) # cacula dQps/dVp dQfdV(graph,0,de)
                if upfc.s in var_v.keys():
                    H[i][var_v[upfc.s]+n_teta]=upfc.dQpsdVs(graph) # cacula dQps/dVs  dPfdV(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  graph[item.k].bar.type!=0:
                    soma2=soma2+upfc.dQspdVs(graph) # calcula dQsp/dVs dfQf(graph,1,para)
                if upfc.p in var_v.keys():
                    H[i][var_v[upfc.p]+n_teta]=upfc.dQspdVp(graph) #cacula dQsp/dVp  dfPfdV(graph,1,de)    
            if  graph[item.k].bar.type!=0 and graph[item.k].bar.type!=1 and (item.k in var_v.keys()):
                if graph[item.k].FlagBS==1:
                    soma2=soma2-2*graph[item.k].Bs*graph[item.k].V ## Colocar as derivadas do shunt## Colocar as derivadas do shunt## Colocar as derivadas do shunt## Colocar as derivadas do shunt
                if graph[item.k].FlagSVC==1:
                    soma2=soma2-2*graph[item.k].SVC.Gk*graph[item.k].SVC.V
                H[i][var_v[item.k]+n_teta]=soma2 # fazer o mesmo para as derivadas do upfc
            soma2=0
        elif item.type==2:
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            #---------------------derivadas em função dos branchs fixos e do TCSC
            if km in graph[k].adjk.keys(): #inserir os ifs dos ramos da upfc
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].adjk[km].dPfdt(graph,0,k)
                if k in var_v.keys():
                    H[i][var_v[k]+n_teta]= graph[k].adjk[km].dPfdV(graph,0,k)
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].adjk[km].dPfdt(graph,0,m)
                if m in var_v.keys():
                    H[i][var_v[m]+n_teta]= graph[k].adjk[km].dPfdV(graph,0,m)
            elif mk in graph[k].adjm.keys(): # entrou como k-m e a barra k é na verdade a m
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].adjm[mk].dPfdt(graph,1,k) # dPmk/dm
                if k in var_v.keys():
                    H[i][var_v[k]+n_teta]= graph[k].adjm[mk].dPfdV(graph,1,k) # dPmk/dm
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].adjm[mk].dPfdt(graph,1,m) # dPmk/dk
                if m in var_v.keys():
                    H[i][var_v[m]+n_teta]= graph[k].adjm[mk].dPfdV(graph,1,m) # dPmk/dk
            #---------------------derivadas em função dos branchs fixos e do UPFC
            elif km in graph[k].bUFPC_adjk.keys(): # medida de fluxo de k-m e o branch entrou como km 
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].bUFPC_adjk[km].dPpsdtp(graph) #dPkmdtk
                if k in var_v.keys():
                    H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjk[km].dPpsVp(graph) #dPkmdtvk
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].bUFPC_adjk[km].dPpsdts(graph) # Dpkmdtm detrivada em relação barra to
                if m in var_v.keys():
                    H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjk[km].dPpsdVs(graph) # 
            elif mk in graph[k].bUFPC_adjm.keys():
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].bUFPC_adjm[mk].dPspdts(graph) 
                if k in var_v.keys():
                    H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjm[mk].dPspdVs(graph)
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].bUFPC_adjm[mk].dPspdtp(graph) # dPpsds
                if m in var_v.keys():
                    H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjm[mk].dPspdVp(graph)
            else:
                print("erro ao calcular fluxo na Jacobiana, medida Fluxo de P {:d}-{:d}".format(graph[k].id,graph[m].id))
                exit(1)
        elif item.type==3:
                k=item.k
                m=item.m
                km=str(k)+"-"+str(m)
                mk=str(m)+"-"+str(k)
                if km in graph[k].adjk.keys(): #inserir os ifs dos ramos da upfc
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].adjk[km].dQfdt(graph,0,k)
                    if k in var_v.keys():    
                        H[i][var_v[k]+n_teta]= graph[k].adjk[km].dQfdV(graph,0,k)
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].adjk[km].dQfdt(graph,0,m)
                    if m in var_v.keys():  
                        H[i][var_v[m]+n_teta]= graph[k].adjk[km].dQfdV(graph,0,m)
                elif mk in graph[k].adjm.keys():
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].adjm[mk].dQfdt(graph,1,k)
                    if k in var_v.keys():
                        H[i][var_v[k]+n_teta]= graph[k].adjm[mk].dQfdV(graph,1,k)
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].adjm[mk].dQfdt(graph,1,m)
                    if m in var_v.keys():
                        H[i][var_v[m]+n_teta]= graph[k].adjm[mk].dQfdV(graph,1,m)
            #---------------------derivadas em função dos branchs fixos e do UPFC
                elif km in graph[k].bUFPC_adjk.keys(): # medida de fluxo de k-m e o branch entrou como km 
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].bUFPC_adjk[km].dQpsdtp(graph) #dPkmdtk
                    if k in var_v.keys():
                        H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjk[km].dQpsVp(graph) #dPkmdtvk
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].bUFPC_adjk[km].dQpsdts(graph) # Dpkmdtm detrivada em relação barra to
                    if m in var_v.keys():
                        H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjk[km].dQpsdVs(graph) # 
                elif mk in graph[k].bUFPC_adjm.keys():
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].bUFPC_adjm[mk].dQspdts(graph) 
                    if k in var_v.keys():
                        H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjm[mk].dQspdVs(graph)
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].bUFPC_adjm[mk].dQspdtp(graph) # dPpsds
                    if m in var_v.keys():
                        H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjm[mk].dQspdVp(graph)
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
        elif item.type==4:
                for key in var_x.keys():
                    H[i][var_x[key]]=0 
        i=i+1

def calc_H_fp_SVC(z,var_svc,graph,H):
    i=0

    for item in z:
        if item.type==0:
            k=item.k
            if graph[k].FlagSVC==1:
                H[i][var_svc[k]]= (graph[k].V**2)*graph[k].SVC.dGkdBsvc()
        if item.type==1:
            k=item.k
            if graph[k].FlagSVC==1:
                H[i][var_svc[k]]= -(graph[k].V**2)*graph[k].SVC.dBkdBsvc()
        i=i+1


def calc_H_fp_UPFC(z,var_UPFC,var_UPFC_vsh,graph,H,Hsh):
    i=0
    n_UPFCs=len(var_UPFC)
    for item in z:
        if item.type==0: #medida de potencia ativa
            k=item.k
            for key, upfc in graph[k].bUFPC_adjk.items():
                H[i][var_UPFC[key]]=upfc.dPpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dPpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dPpsdVse(graph)
                if key in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[key]]=upfc.dPpsdVsh(graph)
            for key, upfc in graph[k].bUFPC_adjm.items():
                H[i][var_UPFC[key]]=upfc.dPspdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dPspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dPspdVse(graph)
                if key in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[key]]=upfc.dPspdVsh(graph)
        if item.type==1: #medida de potencia reativa
            k=item.k
            for key, upfc in graph[k].bUFPC_adjk.items():
                H[i][var_UPFC[key]]=upfc.dQpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dQpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dQpsdVse(graph)
                if key in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[key]]=upfc.dQpsdVsh(graph)
            for key, upfc in graph[k].bUFPC_adjm.items():
                H[i][var_UPFC[key]]=upfc.dQspdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dQspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dQspdVse(graph)
                if key in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[key]]=upfc.dQspdVsh(graph)
        if item.type==2: #medida de potencia ativa
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].bUFPC_adjk.keys():
                H[i][var_UPFC[km]]= graph[k].bUFPC_adjk[km].dPpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dPpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dPpsdVse(graph)
                if km in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[km]]=upfc.dPpsdVsh(graph)
            elif mk in graph[k].bUFPC_adjm.keys():
                H[i][var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dPspdtse(graph)
                H[i][n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dPspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dPspdVse(graph)
                if mk in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[mk]]=upfc.dPspdVsh(graph)
        if item.type==3: #medida de potencia ativa
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].bUFPC_adjk.keys():
                H[i][var_UPFC[km]]= graph[k].bUFPC_adjk[km].dQpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dQpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dQpsdVse(graph)
                if km in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[km]]=upfc.dQpsdVsh(graph)
            elif mk in graph[k].bUFPC_adjm.keys():
                H[i][var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dQspdtse(graph)
                H[i][n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dQspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dQspdVse(graph)
                if mk in var_UPFC_vsh.keys():
                    Hsh[i][var_UPFC_vsh[mk]]=upfc.dQspdVsh(graph)
        i=i+1


def calc_H_EE_UPFC(z,var_UPFC,graph,H):
    i=0
    n_UPFCs=len(var_UPFC)
    for item in z:
        if item.type==0: #medida de potencia ativa
            k=item.k
            for key, upfc in graph[k].bUFPC_adjk.items():
                H[i][var_UPFC[key]]=upfc.dPpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dPpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dPpsdVse(graph)            
                H[i][3*n_UPFCs+var_UPFC[key]]=upfc.dPpsdVsh(graph)
            for key, upfc in graph[k].bUFPC_adjm.items():
                H[i][var_UPFC[key]]=upfc.dPspdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dPspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dPspdVse(graph)
                H[i][3*n_UPFCs+var_UPFC[key]]=upfc.dPspdVsh(graph)
        if item.type==1: #medida de potencia reativa
            k=item.k
            for key, upfc in graph[k].bUFPC_adjk.items():
                H[i][var_UPFC[key]]=upfc.dQpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dQpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dQpsdVse(graph)
                H[i][3*n_UPFCs+var_UPFC[key]]=upfc.dQpsdVsh(graph)
            for key, upfc in graph[k].bUFPC_adjm.items():
                H[i][var_UPFC[key]]=upfc.dQspdtse(graph)
                H[i][n_UPFCs+var_UPFC[key]]=upfc.dQspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[key]]=upfc.dQspdVse(graph)
                H[i][3*n_UPFCs+var_UPFC[key]]=upfc.dQspdVsh(graph)
        if item.type==2: #medida de potencia ativa
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].bUFPC_adjk.keys():
                H[i][var_UPFC[km]]= graph[k].bUFPC_adjk[km].dPpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dPpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dPpsdVse(graph)
                H[i][3*n_UPFCs+var_UPFC[km]]=upfc.dPpsdVsh(graph)
            elif mk in graph[k].bUFPC_adjm.keys():
                H[i][var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dPspdtse(graph)
                H[i][n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dPspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dPspdVse(graph)
                H[i][3*n_UPFCs+var_UPFC[mk]]=upfc.dPspdVsh(graph)
        if item.type==3: #medida de potencia ativa
            k=item.k
            m=item.m
            km=str(k)+"-"+str(m)
            mk=str(m)+"-"+str(k)
            if km in graph[k].bUFPC_adjk.keys():
                H[i][var_UPFC[km]]= graph[k].bUFPC_adjk[km].dQpsdtse(graph)
                H[i][n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dQpsdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[km]]= graph[k].bUFPC_adjk[km].dQpsdVse(graph)
                H[i][3*n_UPFCs+var_UPFC[km]]=upfc.dQpsdVsh(graph)
            elif mk in graph[k].bUFPC_adjm.keys():
                H[i][var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dQspdtse(graph)
                H[i][n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dQspdtsh(graph)
                H[i][2*n_UPFCs+var_UPFC[mk]]= graph[k].bUFPC_adjm[mk].dQspdVse(graph)
                H[i][3*n_UPFCs+var_UPFC[mk]]=upfc.dQspdVsh(graph)
        i=i+1



def calc_C_fp_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,var_UPFC_vsh,graph,C_UPFC):
    i=0

    off=len(var_t)+len(var_v)+len(var_x)+len(var_svc)
    n_upfcs=len(var_UPFC)
    for key in var_UPFC.keys():
        p,s=key.split("-") 
        p=int(p)
        s=int(s)
        C_UPFC[i][var_t[p]]=graph[p].bUFPC_adjk[key].dIgdtp(graph)
        C_UPFC[i][var_t[s]]=graph[p].bUFPC_adjk[key].dIgdts(graph)
        if p in var_v.keys():
            C_UPFC[i][var_v[p]+len(var_t)]=graph[p].bUFPC_adjk[key].dIgdVp(graph)
        if s in var_v.keys():
            C_UPFC[i][var_v[s]+len(var_t)]=graph[p].bUFPC_adjk[key].dIgdVs(graph)
        C_UPFC[i][off+var_UPFC[key]]=graph[p].bUFPC_adjk[key].dIgdtse(graph)
        C_UPFC[i][off+n_upfcs+var_UPFC[key]]=graph[p].bUFPC_adjk[key].dIgdtsh(graph)
        C_UPFC[i][off+2*n_upfcs+var_UPFC[key]]=graph[p].bUFPC_adjk[key].dIgdVse(graph)
        if key in var_UPFC_vsh.keys():
            C_UPFC[i][off+3*n_upfcs+var_UPFC_vsh[key]]=graph[p].bUFPC_adjk[key].dIgdVsh(graph)
        i=i+1


def calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC):
    i=0

    off=len(var_t)+len(var_v)+len(var_x)+len(var_svc)
    n_upfcs=len(var_UPFC)
    for key in var_UPFC.keys():
        p,s=key.split("-") 
        p=int(p)
        s=int(s)
        if p in var_t.keys():
            C_UPFC[i][var_t[p]]=graph[p].bUFPC_adjk[key].dIgdtp(graph)
        if s in var_t.keys():
            C_UPFC[i][var_t[s]]=graph[p].bUFPC_adjk[key].dIgdts(graph)
        C_UPFC[i][var_v[p]+len(var_t)]=graph[p].bUFPC_adjk[key].dIgdVp(graph)
        C_UPFC[i][var_v[s]+len(var_t)]=graph[p].bUFPC_adjk[key].dIgdVs(graph)
        C_UPFC[i][off+var_UPFC[key]]=graph[p].bUFPC_adjk[key].dIgdtse(graph)
        C_UPFC[i][off+n_upfcs+var_UPFC[key]]=graph[p].bUFPC_adjk[key].dIgdtsh(graph)
        C_UPFC[i][off+2*n_upfcs+var_UPFC[key]]=graph[p].bUFPC_adjk[key].dIgdVse(graph)
        C_UPFC[i][off+3*n_upfcs+var_UPFC[key]]=graph[p].bUFPC_adjk[key].dIgdVsh(graph)
        i=i+1



def calc_H_EE_SVC(z,var_svc,graph,H):
    
    i=0
    for item in z:
        if item.type==0:
            k=item.k
            if graph[k].FlagSVC==1:
                H[i][var_svc[k]]= (graph[k].V**2)*graph[k].SVC.dGkdBsvc()
        if item.type==1:
            k=item.k
            if graph[k].FlagSVC==1:
                H[i][var_svc[k]]= -(graph[k].V**2)*graph[k].SVC.dBkdBsvc()
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
             #-------------------ramos de branchs normais e TCSC----------------------------------------#
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
            #-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if  item.k in bar_t:
                    soma1=soma1+upfc.dPpsdtp(graph) # cacula dPps/dtp dPfdt(graph,0,de)
                if upfc.s in bar_t:
                    H[i][var_t[upfc.s]]=upfc.dPpsdts(graph) # cacula dPps/dts  dPfdt(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  item.k in bar_t:
                    soma1=soma1+upfc.dPspdts(graph) # calcula dPsp/dts dfPf(graph,1,para)
                if upfc.p in bar_t:
                    H[i][var_t[upfc.p]]=upfc.dPspdtp(graph) #cacula dPsp/dtp  dfPf(graph,1,de)
            if  item.k in bar_t:
                H[i][var_t[item.k]]=soma1
            soma1=0
            #-------------------ramos de branchs normais e TCSC----------------------------------------#
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
             #-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if item.k in bar_v:
                    soma2=soma2+upfc.dPpsdVp(graph) # cacula dPps/dVp dPfdV(graph,0,de)
                if upfc.s in bar_v:
                    H[i][var_v[upfc.s]+n_teta]=upfc.dPpsdVs(graph) # cacula dPps/dVs  dPfdV(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  item.k in bar_v:
                    soma2=soma2+upfc.dPspdVs(graph) # calcula dPsp/dVs dfPf(graph,1,para)
                if upfc.p in bar_v:
                    H[i][var_v[upfc.p]+n_teta]=upfc.dPspdVp(graph) #cacula dPsp/dVp  dfPfdV(graph,1,de)          
            if item.k in bar_v:
                if graph[item.k].FlagSVC==1:
                    soma2=soma2-2*graph[item.k].SVC.Gk*graph[item.k].V
                if graph[item.k].FlagSVC==1:
                    soma2=soma2+ 2*graph[item.k].SVC.Gk*graph[item.k].V
                H[i][var_v[item.k]+n_teta]=soma2
            soma2=0
        elif item.type==1:
#-------------------ramos de branchs normais e TCSC----------------------------------------#
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
#-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if  item.k in bar_t:
                    soma1=soma1+upfc.dQpsdtp(graph) # cacula dQps/dtp dQfdt(graph,0,de)
                if upfc.s in bar_t:
                    H[i][var_t[upfc.s]]=upfc.dQpsdts(graph) # cacula dQps/dts  dQfdt(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  item.k in bar_t:
                    soma1=soma1+upfc.dQspdts(graph) # calcula dQsp/dts dfQf(graph,1,para)
                if upfc.p in bar_t:
                    H[i][var_t[upfc.p]]=upfc.dQspdtp(graph) #cacula dQsp/dtp  dfPf(graph,1,de)  
            if  item.k in bar_t:
                H[i][var_t[item.k]]=soma1
            soma1=0
            #-------------------ramos de branchs normais e TCSC----------------------------------------#
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
            #-------------------ramos de branchs UPFC----------------------------------------#
            for key,upfc in graph[item.k].bUFPC_adjk.items():# o branch entra com p-s e barra p é a variável
                if  item.k in bar_v:
                    soma2=soma2+upfc.dQpsdVp(graph) # cacula dQps/dVp dQfdV(graph,0,de)
                if upfc.s in bar_v:
                    H[i][var_v[upfc.s]+n_teta]=upfc.dQpsdVs(graph) # cacula dQps/dVs  dPfdV(graph,0,para)
            for key,upfc in graph[item.k].bUFPC_adjm.items(): #o branch entra com p-s e barra s é a variável
                if  item.k in bar_v:
                    soma2=soma2+upfc.dQspdVs(graph) # calcula dQsp/dVs dfQf(graph,1,para)
                if upfc.p in bar_v:
                    H[i][var_v[upfc.p]+n_teta]=upfc.dQspdVp(graph) #cacula dQsp/dVp  dfPfdV(graph,1,de)  
            if graph[item.k].FlagBS==1:
                soma2=soma2-2*graph[item.k].Bs*graph[item.k].V 
            if graph[item.k].FlagSVC==1:
                soma2=soma2-2*graph[item.k].SVC.Bk*graph[item.k].V    
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
            #---------------------derivadas em função dos branchs fixos e do UPFC
            elif km in graph[k].bUFPC_adjk.keys(): # medida de fluxo de k-m e o branch entrou como km 
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].bUFPC_adjk[km].dPpsdtp(graph) #dPkmdtk
                if k in var_v.keys():
                    H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjk[km].dPpsdVp(graph) #dPkmdtvk
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].bUFPC_adjk[km].dPpsdts(graph) # Dpkmdtm detrivada em relação barra to
                if m in var_v.keys():
                    H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjk[km].dPpsdVs(graph) # 
            elif mk in graph[k].bUFPC_adjm.keys():
                if k in var_t.keys():
                    H[i][var_t[k]]= graph[k].bUFPC_adjm[mk].dPspdts(graph) 
                if k in var_v.keys():
                    H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjm[mk].dPspdVs(graph)
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].bUFPC_adjm[mk].dPspdtp(graph) # dPpsds
                if m in var_v.keys():
                    H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjm[mk].dPspdVp(graph)
            
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
                #---------------------derivadas em função dos branchs fixos e do UPFC
                elif km in graph[k].bUFPC_adjk.keys(): # medida de fluxo de k-m e o branch entrou como km 
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].bUFPC_adjk[km].dQpsdtp(graph) #dPkmdtk
                    if k in var_v.keys():
                        H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjk[km].dQpsdVp(graph) #dPkmdtvk
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].bUFPC_adjk[km].dQpsdts(graph) # Dpkmdtm detrivada em relação barra to
                    if m in var_v.keys():
                        H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjk[km].dQpsdVs(graph) # 
                elif mk in graph[k].bUFPC_adjm.keys():
                    if k in var_t.keys():
                        H[i][var_t[k]]= graph[k].bUFPC_adjm[mk].dQspdts(graph) 
                    if k in var_v.keys():
                        H[i][var_v[k]+n_teta]= graph[k].bUFPC_adjm[mk].dQspdVs(graph)
                    if m in var_t.keys():
                        H[i][var_t[m]]= graph[k].bUFPC_adjm[mk].dQspdtp(graph) # dPpsds
                    if m in var_v.keys():
                        H[i][var_v[m]+n_teta]= graph[k].bUFPC_adjm[mk].dQspdVp(graph)
                
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

def calc_cUPFC(graph,var_UPFC,c):
    i=0

    for key in var_UPFC.keys():
        p,s=key.split("-") 
        p=int(p)
        c[i]=-graph[p].bUFPC_adjk[key].Pse(graph)+graph[p].bUFPC_adjk[key].Psh(graph)
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


def new_X_SVC(graph,offset,var_svc,dx):
    for key,item in var_svc.items():
        graph[key].SVC.BSVC=graph[key].SVC.BSVC+dx[offset+item]
        graph[key].SVC.attYk()

def new_X_UPFC(graph,offset,var_UPFC,var_UPFC_sh,dx):
    
    n_upfc=len(var_UPFC)
    for key,item in var_UPFC.items():
        p,s = key.split("-")
        p=int(p)
        graph[p].bUFPC_adjk[key].t_se=graph[p].bUFPC_adjk[key].t_se+dx[offset+var_UPFC[key]]
        graph[p].bUFPC_adjk[key].t_sh=graph[p].bUFPC_adjk[key].t_sh+dx[offset+n_upfc+var_UPFC[key]]
        graph[p].bUFPC_adjk[key].Vse=graph[p].bUFPC_adjk[key].Vse+dx[offset+2*n_upfc+var_UPFC[key]]
        if key in var_UPFC_sh.keys():
            graph[p].bUFPC_adjk[key].Vsh=graph[p].bUFPC_adjk[key].Vsh+dx[offset+3*n_upfc+var_UPFC_sh[key]]

def new_X_EE_UPFC(graph,offset,var_UPFC,dx):
    
    n_upfc=len(var_UPFC)
    for key,item in var_UPFC.items():
        p,s = key.split("-")
        p=int(p)
        graph[p].bUFPC_adjk[key].t_se=graph[p].bUFPC_adjk[key].t_se+dx[offset+var_UPFC[key]]
        graph[p].bUFPC_adjk[key].t_sh=graph[p].bUFPC_adjk[key].t_sh+dx[offset+n_upfc+var_UPFC[key]]
        graph[p].bUFPC_adjk[key].Vse=graph[p].bUFPC_adjk[key].Vse+dx[offset+2*n_upfc+var_UPFC[key]]
        graph[p].bUFPC_adjk[key].Vsh=graph[p].bUFPC_adjk[key].Vsh+dx[offset+3*n_upfc+var_UPFC[key]]





def new_X_TCSC(graph,nvars,var_x,dx):
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        graph[k].adjk[key].xtcsc=graph[k].adjk[key].xtcsc+dx[item+nvars]
        graph[k].adjk[key].AttY()


def checklim_X_TCSC(graph,var_x):
    x_lim_sup=0.2
    x_lim_inf=-0.2
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        if graph[k].adjk[key].xtcsc> x_lim_sup:
            graph[k].adjk[key].xtcsc=x_lim_sup
        elif graph[k].adjk[key].xtcsc< x_lim_inf:
            graph[k].adjk[key].xtcsc=x_lim_inf
        elif np.abs(graph[k].adjk[key].xtcsc)<1e-6:
            if graph[k].adjk[key].xtcsc >0:
                graph[k].adjk[key].xtcsc=1e-6
            else:
                graph[k].adjk[key].xtcsc=-1e-6

        graph[k].adjk[key].AttY()

def new_X_TCSC_lim(graph,nvars,var_x,dx):
    for key,item in var_x.items():
        x_lim_sup_p=0.10
        x_lim_inf_n=-0.10
        x_lim_inf_p=1e-6
        x_lim_sup_n=-1e-6
        k=int(key.split("-")[0])
        # print("{:e}".format(graph[k].adjk[key].xtcsc+dx[item+nvars] ))
        if graph[k].adjk[key].xtcsc+dx[item+nvars] > x_lim_sup_p:
            graph[k].adjk[key].xtcsc= graph[k].adjk[key].xtcsc_ini
        elif graph[k].adjk[key].xtcsc+dx[item+nvars] < x_lim_inf_n:
            graph[k].adjk[key].xtcsc= graph[k].adjk[key].xtcsc_ini
        elif (graph[k].adjk[key].xtcsc+dx[item+nvars] > 0) & (graph[k].adjk[key].xtcsc+dx[item+nvars] < x_lim_inf_p):
            graph[k].adjk[key].xtcsc= x_lim_inf_p
        elif (graph[k].adjk[key].xtcsc+dx[item+nvars] < 0) & (graph[k].adjk[key].xtcsc+dx[item+nvars] > x_lim_sup_n):
            graph[k].adjk[key].xtcsc= x_lim_sup_n
        else:
            graph[k].adjk[key].xtcsc=graph[k].adjk[key].xtcsc+dx[item+nvars]
        # print("xtcsc: {:e}".format(graph[k].adjk[key].xtcsc))
        graph[k].adjk[key].AttY()

def new_X_TCSC_lim2(graph,nvars,var_x,dx):
    for key,item in var_x.items():
        x_lim_inf_p=1e-6
        x_lim_sup_n=-1e-6
        k=int(key.split("-")[0])
        if (graph[k].adjk[key].xtcsc+dx[item+nvars] > 0) & (graph[k].adjk[key].xtcsc+dx[item+nvars] < x_lim_inf_p):
            graph[k].adjk[key].xtcsc= x_lim_inf_p
        elif (graph[k].adjk[key].xtcsc+dx[item+nvars] < 0) & (graph[k].adjk[key].xtcsc+dx[item+nvars] > x_lim_sup_n):
            graph[k].adjk[key].xtcsc= x_lim_sup_n
        else:
            graph[k].adjk[key].xtcsc=graph[k].adjk[key].xtcsc+dx[item+nvars]
        # print("xtcsc: {:e}".format(graph[k].adjk[key].xtcsc))
        graph[k].adjk[key].AttY()


def a_X_TCSC_lim(graph,nvars,var_x,dx):
    x_lim_sup_p=5.000
    x_lim_inf_n=-5.000
    A=[]
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        print("{:e}".format(graph[k].adjk[key].xtcsc+dx[item+nvars] ))
        if graph[k].adjk[key].xtcsc+dx[item+nvars] > x_lim_sup_p:
            a=0.95*(x_lim_sup_p-graph[k].adjk[key].xtcsc)/dx[item+nvars]
        elif graph[k].adjk[key].xtcsc+dx[item+nvars] < x_lim_inf_n:
            a=0.95*(x_lim_inf_n-graph[k].adjk[key].xtcsc)/dx[item+nvars]
        else:
            a=1
        A.append(a)
    return A

def X_TCSC_its(graph,nvars,var_x,dx):

    A=[]
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        dx[item+nvars]=0
def dx_TCSC_max(graph,nvars,var_x,dx):    
    DXs=[]
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        DXs.append(np.max(dx[item+nvars]))
    return np.max(DXs)

def new_X_TCSCC_B(graph,nvars,var_x,dx):
    for key,item in var_x.items():
        k=int(key.split("-")[0])
        graph[k].adjk[key].xtcsc=1/(1/(graph[k].adjk[key].xtcsc)+dx[item+nvars])
        graph[k].adjk[key].AttY()





def load_flow_FACTS(graph,prt=0,tol=1e-12,inici=1,itmax=20):
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
    var_svc=create_x_loadflow_SVC(graph,var_v)
    [z_PUFPC,var_UPFC,var_UPFC_vsh]=create_x_loadflow_UPFC(graph,var_v)


    z=z+zPf+z_PUFPC
    FACTSini(graph,useDFACTS=1)
    Vinici_lf(graph,useDBAR=inici,var_t=var_t,var_x=var_x,z=z)
    # Vinici_DBAR(graph)
    dz=np.zeros(len(z))
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    HUPFC=np.zeros((len(z),3*len(var_UPFC)))
    HUPFC_sh=np.zeros((len(z),len(var_UPFC_vsh)))
    nvar=len(var_v)+len(var_t)+len(var_x)+len(var_svc)+3*len(var_UPFC)+len(var_UPFC_vsh)
    c_UPFC=np.zeros(len(var_UPFC))
    C_UPFC=np.zeros((len(var_UPFC),nvar))

    it=0
    conv=0
    lstdx=[]
    lstdz=[]

    while it<itmax:
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_UPFC)
        calc_H_fp(z,var_t,var_v,graph,H)
        calc_H_fp_TCSC(z,var_x,graph,HTCSC)
        calc_H_fp_SVC(z,var_svc,graph,HSVC)
        calc_H_fp_UPFC(z,var_UPFC,var_UPFC_vsh,graph,HUPFC,HUPFC_sh)
        calc_C_fp_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,var_UPFC_vsh,graph,C_UPFC)

        #calc H fp UPFC
        # calc C fp UPFC

        Hx=np.concatenate((H,HTCSC,HSVC,HUPFC,HUPFC_sh),axis=1)
        Hx=np.concatenate((Hx,C_UPFC),axis=0)
        print("{:e}".format(np.linalg.cond(Hx)))
        A=sparse.csc_matrix(Hx, dtype=float)
        b=np.concatenate((dz,c_UPFC))

       
        dx=sliang.spsolve(A,b)
        if it==0:
            np.savetxt("Hmatrix.csv",Hx,fmt="%.18e",delimiter=",")
            np.savetxt("bvect.csv",b,fmt="%.18e",delimiter=",")
            np.savetxt("dx.csv",dx,fmt="%.18e",delimiter=",")
        
        if 0.1<dx_TCSC_max(graph,len(var_t)+len(var_v),var_x,dx):
            X_TCSC_its(graph,len(var_t)+len(var_v),var_x,dx)    

        new_X(graph,var_t,var_v,dx)
        # if it>2:
        new_X_TCSC_lim(graph,len(var_t)+len(var_v),var_x,dx)
            
        new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,dx)
        new_X_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,var_UPFC_vsh,dx)#
        maxdx=np.max(np.abs(dx))
        maxdz=np.max(np.abs(dz))
        maxc=np.max(np.abs(c_UPFC),initial=0)
        print("max dx {:e} | max dz {:e} | max cUPFC {:e}  ".format(maxdx,maxdz,maxc))
        lstdx.append(maxdx)
        lstdz.append(maxdz)
        if maxdx< tol and maxdz < tol:
            print("convergiu em {} itereacoes".format(it))
            upfc_angle(graph)
            prt_state(graph)
            prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
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
    var_svc=create_x_loadflow_SVC(graph,var_v)

    z=z+zPf
    FACTSini(graph,useDFACTS=1)
    Vinici_lf(graph,useDBAR=inici,var_t=var_t,var_x=var_x,z=z)
    dz=np.zeros(len(z))
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    it=0
    conv=0
    lstdx=[]
    lstdz=[]
    while it<40:
        calc_dz(z,graph,dz)
        calc_H_fp(z,var_t,var_v,graph,H)
        calc_H_fp_TCSC_B(z,var_x,graph,HTCSC)
        calc_H_fp_SVC(z,var_svc,graph,HSVC)
        Hx=np.concatenate((H,HTCSC,HSVC),axis=1)
        A=sparse.csc_matrix(Hx, dtype=float)
        dx=sliang.spsolve(A,dz)
        new_X(graph,var_t,var_v,dx)
        new_X_TCSCC_B(graph,len(var_t)+len(var_v),var_x,dx)
        new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,dx)
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
            if not isinstance(item,meas):
                W[i][i]=1/(prec_virtual**2)
            elif np.abs(item.val)<1e-6:
                W[i][i]=1/(prec_virtual**2)
            else:
                if np.abs((item.sigma))>prec_virtual:
                    W[i][i]=1/(item.sigma**2)
                else:
                    W[i][i]=1/(prec_virtual**2)
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



def upfc_angle(graph):
    for no in graph:
        if no.FlagUPFC==1:
            for key in no.bUFPC_adjk.keys():
                tse=no.bUFPC_adjk[key].t_se
                vse=no.bUFPC_adjk[key].Vse
                j=complex(0,1)
                Vsecom=vse*np.exp(j*tse)
                no.bUFPC_adjk[key].t_se=np.angle(Vsecom)
                no.bUFPC_adjk[key].Vse=np.absolute(Vsecom)
                tsh=no.bUFPC_adjk[key].t_sh
                vsh=no.bUFPC_adjk[key].Vsh
                j=complex(0,1)
                Vshcom=vsh*np.exp(j*tsh)
                no.bUFPC_adjk[key].t_sh=np.angle(Vshcom)
                no.bUFPC_adjk[key].Vsh=np.absolute(Vshcom)