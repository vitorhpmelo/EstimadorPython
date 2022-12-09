import pandas as pd
import numpy as np
from networkcalc import *



def create_dfFluxo(dfDMEDfp,lstFP):
    """
    function to filter the flow measurements in the dfDMEDfp given the branches in the
    lstFP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lstFP: list with the branches with that type of measurement
    @return: dfFLOW: pandas dataframe with the measurements filterd
    """
    dfFLOW=pd.DataFrame()
    for item in lstFP:
        [de,para]=item.split("-")
        dfFLOW=pd.concat([dfFLOW,dfDMEDfp[((dfDMEDfp["type"]==2) |(dfDMEDfp["type"]==3)) &(dfDMEDfp["de"]==int(de)) & (dfDMEDfp["para"]==int(para))]])
    return dfFLOW

def create_dfIP(dfDMEDfp,lst_IP):
    """
    Funcion to filter the power injection measurements in the the dfDMEDfp given buses in the
    list lst_IP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow  and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lst_IP: list with the buses with that type of measurement
    """
    if len(lst_IP)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
    return dfDMEDfp[((dfDMEDfp["type"]==0)|(dfDMEDfp["type"]==1)) & (dfDMEDfp["de"].isin(lst_IP))]

def create_dfV(dfDMEDfp,lst_V):
    """
    Funcion to filter the voltage magnitude measurements in the the dfDMEDfp given buses in the
    list lst_V. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow  and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lst_V: list with the buses with that type of measurement
    """
    return dfDMEDfp[(dfDMEDfp["type"]==4)& (dfDMEDfp["de"].isin(lst_V))]


def create_DMED(sys,prec,graph,ram):
    """
    Creates a DMED file with the measurements according to the "measplan.csv" file
    if, it reads the measurements avaible in the "DMED_fp.csv" file, if it do not exits it runs
    the load flow and creates it. The function recives @sys a string with the name of the system's file
    and the prec dictionary with the pr parameter for each measurement
    @param: sys-string with the name of the system's file
    @param: prec - dictionary with the precision of each measurement type
    @return: dfDMED - pandas dictionary with the measurement set   
    """
    #read the file with the measurement plan
    try: # if the DMED exists the program reads it, this file is not mandatory for power flow 
        dfDMEDfp=pd.read_csv(sys+"/DMED_fp.csv",header=None)
        dfDMEDfp.columns=["type","de","para","zmed","prec"]
    except:
        conv = load_flow(graph,tol=1e-10)
        save_DMED_fp(graph,ram,sys)
        dfDMEDfp=pd.read_csv(sys+"/DMED_fp.csv",header=None)
        dfDMEDfp.columns=["type","de","para","zmed","prec"]
    try:
        df=pd.read_csv(sys+"/measplan.csv",keep_default_na=False)
    except:
        print("There is no measurement plan file")
        exit()

    SCADAlstIP=list(np.int32(list(filter(None,df["PISCADA"].to_list()))))
    SCADAlstFP=list(filter(None,df["PFSCADA"].to_list()))
    SCADAlstV=list(np.int32(list(filter(None,df["VSCADA"].to_list()))))
    SMlstIP=list(np.int32(list(filter(None,df["PISM"].to_list()))))
    SMlstFP=list(filter(None,df["PFSM"].to_list()))
    SMlstV=list(np.int32(list(filter(None,df["VSM"].to_list()))))
    PSEUDOlst=list(np.int32(list(filter(None,df["PSEUDO"].to_list()))))
    Plst=dfDMEDfp[((dfDMEDfp["zmed"]==0.000) & (dfDMEDfp["type"]==0))]["de"].tolist()
    Qlst=dfDMEDfp[((dfDMEDfp["zmed"]==0.000) & (dfDMEDfp["type"]==1))]["de"].tolist()
    Vistuaislst=list(set(Plst).intersection(Qlst))
    Vistuaislst=list(set(Vistuaislst)-set(Vistuaislst).intersection(SCADAlstIP+SMlstIP+PSEUDOlst))

    dfPISCADA=create_dfIP(dfDMEDfp,SCADAlstIP)

    dfPFSCADA=create_dfFluxo(dfDMEDfp,SCADAlstFP)

    dfVSCADA=create_dfV(dfDMEDfp,SCADAlstV)

    dfIPSM=create_dfIP(dfDMEDfp,SMlstIP)

    dfFPSM=create_dfFluxo(dfDMEDfp,SMlstFP)

    dfVSM=create_dfV(dfDMEDfp,SMlstV)

    dfPSEUDO=create_dfIP(dfDMEDfp,PSEUDOlst)

    dfVirtuais=create_dfIP(dfDMEDfp,Vistuaislst)
    dfPFSCADA.loc[:]["prec"]=prec["SCADAPF"]
    dfPISCADA.loc[:]["prec"]=prec["SCADAPI"]
    dfVSCADA.loc[:]["prec"]=prec["SCADAV"]
    dfIPSM.loc[:]["prec"]=prec["SMP"]
    dfFPSM.loc[:]["prec"]=prec["SMP"]
    dfVSM.loc[:]["prec"]=prec["SMV"]
    dfPSEUDO.loc[:]["prec"]=prec["PSEUDO"]
    dfVirtuais.loc[:]["prec"]=prec["VIRTUAL"]

    dfDMED=pd.concat([dfPISCADA,dfIPSM,dfPSEUDO,dfVirtuais,dfPFSCADA,dfFPSM,dfVSCADA,dfVSM])
    return dfDMED

def insert_res(dfDMEDsr):
    """
    Inserts gaussian noise in the measurement set, with variance according with the 
    precision and the magnitude of the measurement.
    """
    e=np.random.normal(size=(len(dfDMEDsr)))
    for i in range(len(e)):
        if e[i]>3:
            e[i]=3
        elif e[i]<-3:
            e[i]=-3
    dfDMEDr=dfDMEDsr.copy()
    dfDMEDr.iloc[:]["zmed"]=dfDMEDsr["zmed"]+e*dfDMEDsr["prec"]*np.abs(dfDMEDsr["zmed"])/3
    return dfDMEDr