from csv import DictReader
from classes import *
import pandas as pd
import numpy as np

def creat_bar(dfDBAR):  
    """
    Function to read the information about the network in the data frame and put it into the bar struture
    @param: dfDBAR - data frame with the informations about the network
    @return: vector of instances bar class, with the information about the netowrk buses
    @return: i - number of buses
    @return: pv - index of the pv buses
    @return: pq - index of the pq buses
    @return: ind_id - dict to convert the name of the bus to its index
    """
    bars=[]
    i=0
    nvar=0
    pv=[]
    pq=[]
    ind_id={}
    for idx, row in dfDBAR.iterrows():#reads each line of the data frame
        item=bar(int(row["id"]),int(row["type"]),i) 
        item.V=row.V
        item.teta=row.teta*np.pi/180 # converts the angle to radians
        item.Pg=row["Pg"]/100#converts the power from MW to p.u.
        item.Qg=row["Qg"]/100#converts the power from MW to p.u.
        item.Pd=row["Pd"]/100#converts the power from MW to p.u.
        item.Qd=row["Qd"]/100#converts the power from MW to p.u.
        item.Bs=row["Bs"]/100#converts the power from MW to p.u.
        if int(row["type"])==1:
            pv.append(i)#create the list of PV buses
        elif int(row["type"])==2:
            pq.append(i)#create the list of PQ buses
        bars.append(item)
        ind_id[int(row["id"])]=i # create the dictionary relating the name of the bus to it index
        i=i+1
    return bars,i,pv,pq,ind_id

def create_bran(dfDBRAN,ind_i):
    """
    Function to read the information in the data frame dfDBRAN and put it in the ram dictionary, that is composed by instances of the
    branch class. This ditc contains all the information about the network branches.
    @param: dfDBRAN - Data Frame with the information about the buses
    @param: ind_i - dict to translate the name of the bus to it index
    @return ram - list of instances of branch class with the information about the network branches
    @return i - number of branches 
    """
    ram={}
    i=0
    d={}
  

  ##determina linhas paralelas
    dparalelas={}
    for idx, row in dfDBRAN.iterrows():
        mask1=(dfDBRAN["de"]== row["de"]) & (dfDBRAN["para"]== row["para"])
        mask2=(dfDBRAN["de"]== row["para"]) & (dfDBRAN["para"]== row["de"])
        if sum(mask1) + sum(mask2)>1:
            if str(row["de"])+"-"+str(row["para"]) in dparalelas.keys():
                dparalelas[str(row["de"])+"-"+str(row["para"])].append(idx)
            elif str(row["para"])+"-"+str(row["de"]) in dparalelas.keys():
                dparalelas[str(row["para"])+"-"+str(row["de"])].append(idx)
            else:
                dparalelas[str(row["de"])+"-"+str(row["para"])]=[idx]

    # remove linhas paralelas
    for key,item in dparalelas.items():
        bsh=0
        y=0
        for line in item:
            y=y+1/complex(dfDBRAN.loc[line,"r"],dfDBRAN.loc[line,"x"])
            bsh=bsh+dfDBRAN.loc[line,"bsh"]
        dfDBRAN.loc[item[0],"r"]=np.real(1/y)
        dfDBRAN.loc[item[0],"x"]=np.imag(1/y)
        dfDBRAN.loc[item[0],"bsh"]=bsh

        dfDBRAN.drop(item[1:],inplace=True)
        
    dfDBRAN.reindex()

    for id, row in dfDBRAN.iterrows():
        key=str(ind_i[int(row["de"])])+"-"+str(ind_i[int(row["para"])])
        item=branch(int(row["id"]),ind_i[int(row["de"])],ind_i[int(row["para"])],int(row["type"]),i)
        item.x=row["x"]
        item.r=row["r"]
        item.bsh=complex(0,row["bsh"]/2) #divides the shunt suceptance by two
        item.tap=row["tap"]
        item.cykm()#calculates the ykm
        item.twoPortCircuit()#creates the two port circuit
        ram[key]=item    
        i+=1
    return ram,i



def create_TCSC(dfFACTS,ind_i):
    """
    Function to read the information in the data frame dfTCSCS and put it in the ram dictionary, that is composed by instances of the
    branch class. This ditc contains all the information about the network branches.
    @param: dfFACTS - Data Frame with the information about the buses
    @param: ind_i - dict to translate the name of the bus to it index
    @return ram - list of instances of branch class with the information about the network branches
    @return i - number of branches 
    """
    ram={}
    i=0
    d={}
##determina linhas paralelas
    if dfFACTS.empty:
        ram=[]
        i=0
        return ram,i    
    if dfFACTS[dfFACTS["type"]==0].empty:
        ram=[]
        i=0
        return ram,i
    dfTCSC=dfFACTS[dfFACTS["type"]==0].copy()
    for id, row in dfTCSC.iterrows():
        #ram type 3 == TCSC
        key=str(ind_i[int(row["de"])])+"-"+str(ind_i[int(row["para"])])
        item=branTCSC(int(row["id"]),ind_i[int(row["de"])],ind_i[int(row["para"])],3,i,row["a"],row["xtscc_ini"],row["Pfesp"])
        ram[key]=item    
        i+=1
    return ram,i


def create_SVC(dfFACTS,ind_i):
    """
    Function to read the information in the data frame dfSVC and put it in the SVC dictionary, that is composed by instances of the
    SVC class. This ditc contains all the information about the network SVCs.
    @param: dfFACTS - Data Frame with the information about the buses
    @param: ind_i - dict to translate the name of the bus to it index
    @return ram - list of instances of branch class with the information about the network branches
    @return i - number of branches 
    """

    if dfFACTS.empty:
        svc=[]
        i=0
        return svc,i    

    if dfFACTS[dfFACTS["type"]==1].empty:
        svc=[]
        i=0
        return svc,i
    
    svc={}
    i=0
    d={}
    dfSVC=dfFACTS[dfFACTS["type"]==1].copy()
    for id, row in dfSVC.iterrows():
        #ram type 3 == TCSC
        key=ind_i[int(row["de"])]
        item=SVC(int(row["id"]),ind_i[int(row["de"])],row["Rt"],row["Xt"],row["Bini"],row["Bmax"],row["Bmin"],row["aini"],row["amax"],row["amin"])
        svc[key]=item    
        i+=1
    return svc,i


def create_UPFC(dfFACTS,ind_i):
    """
    Function to read the information in the data frame dfFACTS and exatract the ones about the UPFC and put it into a dictionary, that is composed by instances of the
    UPFC class. This ditc contains all the information about the network UFPCS.
    @param: dfFACTS - Data Frame with the information about the buses
    @param: ind_i - dict to translate the name of the bus to it index
    @return ramUPFC - list of instances of branch class with the information about the network branches
    @return i - number of UPFC 
    """
    ram={}
    i=0
    d={}
##determina linhas paralelas
    if dfFACTS.empty:
        ram={}
        i=0
        return ram,i    
    if dfFACTS[dfFACTS["type"]==2].empty:
        ram={}
        i=0
        return ram,i
    
    dfUPFC=dfFACTS[dfFACTS["type"]==2].copy()
    for id, row in dfUPFC.iterrows():
        #ram type 3 == TCSC

        key=str(ind_i[int(row["from"])])+"-"+str(ind_i[int(row["to"])])
        item=UPFC(id=int(row["id"]),de=ind_i[int(row["from"])],para=ind_i[int(row["to"])],Vse_ini=row["Vse"],\
                t_se_ini=row["t_se"],Vsh_ini=row["Vsh"],t_sh_ini=row["t_sh"],Psp=row["Psp"],Qsp=row["Qsp"],\
                Vp=row["Vp"],Rse=row["Rse"],Xse=row["Xse"],Rsh=row["Rsh"],Xsh=row["Xsh"],Vse_max=row["Vse_max"],\
                Vse_min=row["Vse_min"],Vsh_max=row["Vsh_max"],Vsh_min=row["Vsh_min"],mode=row["mode"])
        ram[key]=item    
        i+=1
    return ram,i

def create_graph(bars,ram):
    """
    Creates the network graph usinf the information of the bars list and the ram dic.
    @param bars - list of instances of the class bar with the information about the network buses
    @param ram - list of instances of branch class with the information about the network branches
    @return graph - list of instances of the class node, that form the network graph
    """
    graph=[]
    for item in bars:
        node=node_graph(item.i,item)
        if abs(item.Bs)>1e-8:
            node.FlagBS=1#informs if exists a capacitor bank or a reactor
            node.Bs=item.Bs
        graph.append(node)
    for key,item in ram.items(): #save the adjacent buses in the node and the rams connected to it
        k=int(key.split("-")[0])
        m=int(key.split("-")[1])
        graph[k].adjk.update({key:item})
        graph[k].ladjk.append(m)
        graph[m].adjm.update({key:item})
        graph[m].ladjm.append(k)
    
    return graph

def addTCSCingraph(graph,ramfacts):

    if not ramfacts:
        return

    for key,item in ramfacts.items(): #save the adjacent buses in the node and the rams connected to it
        if item.type==3:#Type 3 == TCSC
            k=int(key.split("-")[0]) #bus from
            m=int(key.split("-")[1]) #bus to
            item.AttY() # creates adimitance matrix
            graph[k].FlagTCSC=1 # indicates that there is TCSC in the bus from
            graph[m].FlagTCSC=1# indicates that there is TCSC connected in the bus to
            graph[k].bFACTS_adjk.update({key:item}) #inserts the ram key in the bus adj of the bus from dic only of FACTS
            graph[m].bFACTS_adjm.update({key:item}) #inserts the ram key in the bus adj of the bus to dic only of FACTS
            graph[k].adjk.update({key:item}) #inserts the ram key in the bus from adjk dict
            graph[k].ladjk.append(m) #inserts the ram key in the bus from adjk list
            graph[m].adjm.update({key:item}) #inserts the ram key in the bus to adjm dict
            graph[m].ladjm.append(k) #inserts the ram key in the bus to adjm list

def addSVCingraph(graph,busSVC):
     
    if not busSVC:
        return

    for k,svc in busSVC.items():
        graph[k].SVC=svc # it is not a reference of the original object in busSVC
        graph[k].FlagSVC=True


def addUPFCingraph(graph,ramUPFC):

    if not ramUPFC:
        return

    for key,item in ramUPFC.items(): #save the adjacent buses in the node and the rams connected to it
        k=int(key.split("-")[0]) #bus from
        m=int(key.split("-")[1]) #bus to
        graph[k].FlagUPFC=1 # indicates that there is TCSC in the bus from
        graph[m].FlagUPFC=1# indicates that there is TCSC connected in the bus to
        graph[k].bUFPC_adjk.update({key:item}) #inserts the ram key in the bus adj of the bus from dic only of FACTS
        graph[m].bUFPC_adjm.update({key:item}) #inserts the ram key in the bus adj of the bus to dic only of FACTS
