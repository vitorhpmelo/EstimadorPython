from classes import *
import pandas as pd
import numpy as np

def read_files(sys):
    """
    This function performs the reading of the files with the information about the Power System.
    @param: sys - string with the name of the system
    @return dfDBAR - Data Frame with the information about the system's buses
    @return dfDBRAN - Data Frame with the information about the system's branches
    @return dfDMED - Data Frame with the information about the system measurements
    """
    try: # if the DBAR exists the program reads it, if not it stops. This file is mandatory 
        dfDBAR=pd.read_csv(sys+"/DBAR.csv",header=None,dtype={0:np.int64,1:np.int64})
        dfDBAR.columns=["id","type","V","teta","Pg","Qg","Pd","Qd","Bs"]
    except:
        print("Error while reading DBAR file")
        quit()

    try: # if the DBRAN exists the program reads it, if not it stops. This file is mandatory
        dfDBRAN=pd.read_csv(sys+"/DBRAN.csv",header=None,dtype={0:np.int64,1:np.int64,2:np.int64,3:np.int64})
        dfDBRAN.columns=["id","type","de","para","r","x","bsh","tap"]
    except:
        print("Error while reading DBRAN file")
        exit(1)
    try: # if the DMED exists the program reads it, this file is not mandatory for power flow 
        dfDMED=pd.read_csv(sys+"/DMED.csv",header=None)
        dfDMED.columns=["type","de","para","zmed","prec"]
    except:
        print("There is no DMED")
        dfDMED=[]
    try: # if the DMED exists the program reads it, this file is not mandatory for power flow 
        dfDTCSC=pd.read_csv(sys+"/DTCSC.csv",header=None)
        dfDTCSC.columns=["id","de","para","a","xtscc_ini","Pfesp"]
        dfDTCSC["type"]=0
    except:
        print("There is no DTCSC")
        dfDTCSC=[]    
    try:
        dfDFACTS=pd.concat([dfDTCSC], axis=0, ignore_index=True)
    except:
        dfDFACTS=[]
    return dfDBAR,dfDBRAN,dfDMED,dfDFACTS


def prt_state(graph):
    """
    Function to print in the scream the value of the state variables, in the network's graph
    @param: graph Graph structure with the information about the network
    """
    for no in graph:
        s="Barra: {:d} | V : {:f} | teta : {:f}".format(no.bar.id,no.V,no.teta*180/np.pi)
        print(s)


def save_DMED_fp(graph,ram,sys):
    """
    Function to save the file with all measurements possible, from a load flow simulation. 
    It use the graph of the network to calculate every possible measurement and save it in a file called 
    DMED_fp.csv into the system's folder.
    @param: graph Graph structure with the information about the network
    @param: ram dictionary with the information about the network branches
    @sys: string with the system folder's name
    """
    medidas=[] 
    Pinj=[]
    Qinj=[]
    Vmod=[]
    Pkm=[]
    Pmk=[]
    Qkm=[]
    Qmk=[]
    #calculates the Power Inejection (Reactive and Active)
    for no in graph:
        linha=[0,no.bar.id,-1,no.P(graph),1]
        Pinj.append(linha)
        linha=[1,no.bar.id,-1,no.Q(graph),1]
        Qinj.append(linha)
        linha=[4,no.bar.id,-1,no.V,1]
        Vmod.append(linha)

    #calculates the flows in the branches
    for key,r in ram.items():
        #calculate from k to m
        linha=[2,graph[r.de].bar.id,graph[r.para].bar.id,r.Pf(graph,0),1.0]
        linha2=[3,graph[r.de].bar.id,graph[r.para].bar.id,r.Qf(graph,0),1.0]
        Pkm.append(linha)
        Qkm.append(linha2)
        #calculate from m to k
        linha=[2,graph[r.para].bar.id,graph[r.de].bar.id,r.Pf(graph,1),1.0]
        linha2=[3,graph[r.para].bar.id,graph[r.de].bar.id,r.Qf(graph,1),1.0]
        Pmk.append(linha)
        Qmk.append(linha2)

    medidas=Pinj+Qinj+Pkm+Qkm+Pmk+Qmk+Vmod
    dfDMED=pd.DataFrame(medidas,columns=["type","de","para","zmed","pre"])
    dfDMED.to_csv(sys+"/DMED_fp.csv",index=False,float_format="%.7f",header=False)



def save_DBAR(graph):

    id=[]
    tipo=[]
    V=[]
    teta=[]
    Pg=[]
    Qg=[]
    Pd=[]
    Qd=[]
    Bs=[]
    for no in graph:
        id.append(no.bar.id)
        tipo.append(no.bar.type)
        V.append(no.V)
        teta.append(no.teta*180/np.pi)
        Pg.append(no.bar.Pg*100)
        Qg.append(no.bar.Qg*100)
        Pd.append(no.bar.Pd*100)
        Qd.append(no.bar.Qd*100)
        Bs.append(no.bar.Bs*100)


    d={"id":id,"tipo":tipo,"V":V,"teta":teta,"Pg":Pg,"Qg":Qg,"Pd":Pd,"Qd":Qd,"Bs":Bs}
    dfDBAR=pd.DataFrame(d)

    dfDBAR.to_csv("DBAR.csv",header=None,index=None,float_format="%.7f")


