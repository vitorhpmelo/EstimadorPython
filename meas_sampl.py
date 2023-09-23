import pandas as pd
import numpy as np
from networkcalc import *


def create_dfmeasTCSC(dfDMEDfp,lstTCSC):
    """
    function to filter the TCSC measurements in the dfDMEDfp given the branches in the
    lstFP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lstFP: list with the branches with that type of measurement
    @return: dfFLOW: pandas dataframe with the measurements filterd
    """
    dfTCSC=pd.DataFrame()

    if len(lstTCSC)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
    for item in lstTCSC:
        [de,para]=item.split("-")
        dfTCSC=pd.concat([dfTCSC,dfDMEDfp[((dfDMEDfp["type"]==10)) &(dfDMEDfp["de"]==int(de)) & (dfDMEDfp["para"]==int(para))]])
    return dfTCSC



def create_dfmeasUPFCVsh(dfDMEDfp,lstUPFC):
    """
    function to filter the UPFC Vsh measurements in the dfDMEDfp given the branches in the
    lstFP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lstUPFC: list with the branches with that type of measurement
    @return: dfmeasUPFCVsh: pandas dataframe with the measurements filterd
    """
    dfmeasUPFCVsh=pd.DataFrame()

    if len(lstUPFC)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
    for item in lstUPFC:
        [de,para]=item.split("-")
        dfmeasUPFCVsh=pd.concat([dfmeasUPFCVsh,dfDMEDfp[((dfDMEDfp["type"]==12)) &(dfDMEDfp["de"]==int(de)) & (dfDMEDfp["para"]==int(para))]])
    return dfmeasUPFCVsh



def create_dfmeasUPFCtsh(dfDMEDfp,lstUPFC):
    """
    function to filter the UPFC tsh measurements in the dfDMEDfp given the branches in the
    lstFP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lstUPFC: list with the branches with that type of measurement
    @return: dfmeasUPFCtsh: pandas dataframe with the measurements filterd
    """
    dfmeasUPFCtsh=pd.DataFrame()

    if len(lstUPFC)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
    for item in lstUPFC:
        [de,para]=item.split("-")
        dfmeasUPFCtsh=pd.concat([dfmeasUPFCtsh,dfDMEDfp[((dfDMEDfp["type"]==13)) &(dfDMEDfp["de"]==int(de)) & (dfDMEDfp["para"]==int(para))]])
    return dfmeasUPFCtsh

def create_dfmeasUPFCVse(dfDMEDfp,lstUPFC):
    """
    function to filter the UPFC Vse measurements in the dfDMEDfp given the branches in the
    lstFP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lstUPFC: list with the branches with that type of measurement
    @return: dfmeasUPFCVse: pandas dataframe with the measurements filterd
    """
    dfmeasUPFCVse=pd.DataFrame()

    if len(lstUPFC)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
    for item in lstUPFC:
        [de,para]=item.split("-")
        dfmeasUPFCVse=pd.concat([dfmeasUPFCVse,dfDMEDfp[((dfDMEDfp["type"]==14)) &(dfDMEDfp["de"]==int(de)) & (dfDMEDfp["para"]==int(para))]])
    return dfmeasUPFCVse

def create_dfmeasUPFCtse(dfDMEDfp,lstUPFC):
    """
    function to filter the UPFC tse measurements in the dfDMEDfp given the branches in the
    lstFP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lstUPFC: list with the branches with that type of measurement
    @return: create_dfmeasUPFCtse: pandas dataframe with the measurements filterd
    """
    dfmeasUPFCVse=pd.DataFrame()

    if len(lstUPFC)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
    for item in lstUPFC:
        [de,para]=item.split("-")
        dfmeasUPFCVse=pd.concat([dfmeasUPFCVse,dfDMEDfp[((dfDMEDfp["type"]==15)) &(dfDMEDfp["de"]==int(de)) & (dfDMEDfp["para"]==int(para))]])
    return dfmeasUPFCVse



def create_dfmeasSVC(dfDMEDfp,lst_svc):
    """
    Funcion to filter the SVC variable measurements in the the dfDMEDfp given buses in the
    list lst_IP. It recives the dfDMEDfp data frame with all the possible measurements avaible in the network
    obatined by the load flow  and returns only the ones desired.
    @param: dfDMEDfp: pandas dataframe with all the measurements avaible in the loadflow
    @param: lst_svc: list with the buses with that type of measurement
    """
    if len(lst_svc)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
    return dfDMEDfp[((dfDMEDfp["type"]==11)) & (dfDMEDfp["de"].isin(lst_svc))]



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

    if len(lstFP)<1:
        return dfDMEDfp[dfDMEDfp['type']==-1]
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


def create_DMED(sys,prec,graph,ram,dUPFC={},dfDMEDfp=pd.DataFrame()):
    """
    Creates a DMED file with the measurements according to the "measplan.csv" file
    if, it reads the measurements avaible in the "DMED_fp.csv" file, if it do not exits it runs
    the load flow and creates it. The function recives @sys a string with the name of the system's file
    and the prec dictionary with the pr parameter for each measurement
    @param: sys-string with the name of the system's file
    @param: prec - dictionary with the precision of each measurement type
    @return: dfDMED - pandas dictionary with the measurement set   
    """
    #read the file with the measurement pla
    if dfDMEDfp.empty:

        try: # if the DMED exists the program reads it, this file is not mandatory for power flow 
            dfDMEDfp=pd.read_csv(sys+"/DMED_fp.csv",header=None)
            dfDMEDfp.columns=["type","de","para","zmed","prec"]
        except:
            conv = load_flow(graph,tol=1e-10)
            save_DMED_fp(graph,ram,sys,dUPFC)
            dfDMEDfp=pd.read_csv(sys+"/DMED_fp.csv",header=None)
            dfDMEDfp.columns=["type","de","para","zmed","prec"]

    try:
        df=pd.read_csv(sys+"/measplan.csv",keep_default_na=False)
    except:
        print("There is no measurement plan file")
        exit()
    
    prec_standard={"SCADAPF":0.01,"SCADAPI":0.01,"SCADAV":0.01,"SMP":0.01,"SMP":0.01,"SMV":0.01,"PSEUDO":0.01,"VIRTUAL":0.01}
    


    for key in list(set(prec_standard.keys())-(prec.keys())):
        prec[key]=prec_standard[key]

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
    dfPFSCADA.loc[:,"prec"]=prec["SCADAPF"]
    dfPISCADA.loc[:,"prec"]=prec["SCADAPI"]
    dfVSCADA.loc[:,"prec"]=prec["SCADAV"]
    dfIPSM.loc[:,"prec"]=prec["SMP"]
    dfFPSM.loc[:,"prec"]=prec["SMP"]
    dfVSM.loc[:,"prec"]=prec["SMV"]
    dfPSEUDO.loc[:,"prec"]=prec["PSEUDO"]
    dfVirtuais.loc[:,"prec"]=prec["VIRTUAL"]

    dfDMED=pd.concat([dfPISCADA,dfIPSM,dfPSEUDO,dfVirtuais,dfPFSCADA,dfFPSM,dfVSCADA,dfVSM])
    return dfDMED

def insert_res(dfDMEDsr,N=100):
    """
    Inserts gaussian noise in the measurement set, with variance according with the 
    precision and the magnitude of the measurement.
    """
    # np.random.seed(N)
    e=np.random.normal(size=(len(dfDMEDsr)))
    for i in range(len(e)):
        if e[i]>3:
            e[i]=2.0
        elif e[i]<-3:
            e[i]=-3.0
    dfDMEDr=dfDMEDsr.copy()
    dfDMEDr.loc[:,"zmed"]=dfDMEDsr["zmed"]+e*dfDMEDsr["prec"]*np.abs(dfDMEDsr["zmed"])/3
    return dfDMEDr




def create_DMED_FACTS(sys,prec,graph,ram,ramUPFC,dfDMEDfp=pd.DataFrame()):
    """
    Creates a DMED part for the FACTS with the measurements according to the "measplan_FACTS.csv" file
    if, it reads the measurements avaible in the "DMED_fp.csv" file, if it do not exits it runs
    the load flow and creates it. The function recives @sys a string with the name of the system's file
    and the prec dictionary with the pr parameter for each measurement
    @param: sys-string with the name of the system's file
    @param: prec - dictionary with the precision of each measurement type
    @return: dfDMED - pandas dictionary with the measurement set   
    """
    #read the file with the measurement plan
    if dfDMEDfp.empty:
        try: # if the DMED exists the program reads it, this file is not mandatory for power flow 
            dfDMEDfp=pd.read_csv(sys+"/DMED_fp.csv",header=None)
            dfDMEDfp.columns=["type","de","para","zmed","prec"]
        except:
            conv = load_flow(graph,tol=1e-10)
            save_DMED_fp(graph,ram,sys,ramUPFC)
            dfDMEDfp=pd.read_csv(sys+"/DMED_fp.csv",header=None)
            dfDMEDfp.columns=["type","de","para","zmed","prec"]
    try:
        df=pd.read_csv(sys+"/measplanFACTS.csv",keep_default_na=False)
    except:
        print("There is no measurement fatcs plan file")
        quit()


    
    prec_standard={"TCSCvar":0.01,"SVCvar":0.01,"UPFCt_sh":0.01,"UPFCV_sh":0.01,"UPFCt_se":0.01,"UPFCV_se":0.01}
    
    for key in list(set(prec_standard.keys())-(prec.keys())):
        prec[key]=prec_standard[key]


    TCSCvar=list(filter(None,df["TCSCvar"].to_list()))
    SVCvar=list(np.int32(list(filter(None,df["SVCvar"].to_list()))))
    UPFCt_se=list(filter(None,df["UPFCt_se"].to_list()))
    UPFCV_se=list(filter(None,df["UPFCV_se"].to_list()))
    UPFCt_sh=list(filter(None,df["UPFCt_sh"].to_list()))
    UPFCV_sh=list(filter(None,df["UPFCV_sh"].to_list()))


    dfmeasTCSC=create_dfmeasTCSC(dfDMEDfp,TCSCvar)

    dfmeasSVCvar=create_dfmeasSVC(dfDMEDfp,SVCvar)

    dfcreate_dfmeasUPFCVsh=create_dfmeasUPFCVsh(dfDMEDfp,UPFCV_sh)
    dfcreate_dfmeasUPFCtsh=create_dfmeasUPFCtsh(dfDMEDfp,UPFCt_sh)
    dfcreate_dfmeasUPFCVse=create_dfmeasUPFCVse(dfDMEDfp,UPFCV_se)
    dfcreate_dfmeasUPFCtse=create_dfmeasUPFCtse(dfDMEDfp,UPFCt_se)


    dfmeasTCSC.loc[:,"prec"]=prec["TCSCvar"]
    dfmeasSVCvar.loc[:,"prec"]=prec["SVCvar"]
    dfcreate_dfmeasUPFCVsh.loc[:,"prec"]=prec["UPFCt_sh"]
    dfcreate_dfmeasUPFCtsh.loc[:,"prec"]=prec["UPFCV_sh"]
    dfcreate_dfmeasUPFCVse.loc[:,"prec"]=prec["UPFCt_se"]
    dfcreate_dfmeasUPFCtse.loc[:,"prec"]=prec["UPFCV_se"]
 

    dfDMED=pd.concat([dfmeasTCSC,dfmeasSVCvar,dfcreate_dfmeasUPFCVsh,dfcreate_dfmeasUPFCtsh,dfcreate_dfmeasUPFCVse,dfcreate_dfmeasUPFCtse])
    return dfDMED