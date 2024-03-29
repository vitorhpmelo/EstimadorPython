from classes import *
import numpy as np
import pandas as pd
from readfiles import *
import scipy.sparse.linalg as sliang 
import scipy.sparse as sparse 

def Vinici(graph,flatStart=0):
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
    if flatStart==1:
        for no in graph:
            no.V=1
            no.teta=graph[idxref].bar.teta
    else:
        for no in graph:
            no.V=1
            no.teta=0

def Vinici_lf(graph,useDBAR=1):
    '''
    Function to initate the voltages (state variables) for the load flow, 
    PQ buses recive 1 for the voltage module and 0 for the angle,
    PV recive the V from the DBAR for the module
    slack initate with the voltage from the DB 
    @param: graph list of instances of the node class with all the information about the network
    '''
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
                    no.teta=0  
        else:
            no.V=1
            no.teta=0



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
                H[i][var_v[k]+n_teta]= graph[k].adjk[km].dPfdV(graph,0,k)
                if m in var_t.keys():
                    H[i][var_t[m]]= graph[k].adjk[km].dPfdt(graph,0,m)
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


def load_flow(graph,prt=0,tol=1e-6):
    Vinici_lf(graph)
    #Vinici(graph,flatStart=0)
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