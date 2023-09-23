#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from classes import *
from readfiles import *
from networkstruc import *
from SS import *
from meas_sampl import *
import pandas as pd
import numpy as np
from networkcalc import *
from BadData import *
import numpy.linalg as liang
import scipy.sparse.linalg as sliang 
import copy


#%% Lê arquivos e constroi a estrutura da rede

sys="IEEE118_rakp2009"
measFACTS=False

if measFACTS==True:
    Meas="ComMedidas"
else:
    Meas="SemMedidas"



dfDBAR,dfDBRAN,dfDMED,dfDFACTS=read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

[ramTCSC,nbranTCSC]=create_TCSC(dfDFACTS,ind_i)

[busSVC,BUS_SVC]=create_SVC(dfDFACTS,ind_i)

[ramUPFC,nbranUPFC]=create_UPFC(dfDFACTS,ind_i)

graph=create_graph(bars,ram)

addTCSCingraph(graph,ramTCSC)

addSVCingraph(graph,busSVC)

addUPFCingraph(graph,ramUPFC)


#%%

dfTCSC_original_values={}
dfsvc_original_values={}
dfUPFC_original_values={}
for key ,tcsc in ramTCSC.items():
    dfTCSC_original_values[key]=tcsc.Pfesp
for key,svc in busSVC.items():
    dfsvc_original_values[key]=graph[key].bar.V
for key,upfc in ramUPFC.items():
    dfUPFC_original_values[key]={}
    dfUPFC_original_values[key]["Psp"]=upfc.Psp_set
    dfUPFC_original_values[key]["Qsp"]=upfc.Qsp_set
    dfUPFC_original_values[key]["Vp"]=graph[upfc.p].bar.V
        

dfcasos=pd.DataFrame(data={"TCSC":[-15,-10,10,15],"SVC":[-1,-1,1,1],"UPFC_flow":[15,10,-10,-15],"UPFC_V":[-2,-1,1,1],"TCSC_ini":[0.01,0.01,-0.01,-0.01],"SVC_ini":[0.10,0.10,-1.0,-1.0]})
#%%
dDMEDfps={}
dState_ref={}
dStateTCSC_ref={}
for idx, row in dfcasos.iterrows():
    for key ,tcsc in ramTCSC.items():
        tcsc.Pfesp=dfTCSC_original_values[key]*(1+(row["TCSC"]/100))
        tcsc.xtcsc_ini=row["TCSC_ini"]
    for key,svc in busSVC.items():
        graph[key].bar.V=dfsvc_original_values[key]*(1+(row["SVC"]/100))
        svc.Bini=row["SVC_ini"]
    for key,upfc in ramUPFC.items():
        upfc.Psp_set=dfUPFC_original_values[key]["Psp"]*(1+(row["UPFC_flow"]/100))
        upfc.Qsp_set=dfUPFC_original_values[key]["Psp"]*(1+(row["UPFC_flow"]/100))
        graph[upfc.p].bar.V=dfUPFC_original_values[key]["Vp"]*(1+(row["UPFC_V"]/100))

    try:    
        conv=load_flow_FACTS(graph,inici=1,prt=1,itmax=20,printgrad=0,printres=0)
    except:
        conv=0

    #get states and 
    if conv==1:
        ram.update(ramTCSC)
        dDMEDfps[idx]=save_DMED_fp(graph,ram,sys,ramUPFC)
        dState_ref[idx]=get_state(graph)
        dStateTCSC_ref[idx]=get_state_FACTS(ramTCSC,busSVC,ramUPFC)

#%%
dfSATES_ref=pd.DataFrame()
dfSATES_FACTS_ref=pd.DataFrame()

for key,item in dState_ref.items():
    df=item
    df["method"]="PF"
    df["scenario"]=key
    dfSATES_ref=pd.concat([dfSATES_ref,df])
    
for key,item in dStateTCSC_ref.items():
    df=item
    df["method"]="PF"
    df["scenario"]=key
    dfSATES_FACTS_ref=pd.concat([dfSATES_FACTS_ref,df])

dfSATES_ref.to_csv("ResultadosPSCC/state_ref"+sys+".csv")
dfSATES_FACTS_ref.to_csv("ResultadosPSCC/state_FACTS_ref"+sys+".csv")

#%%
if measFACTS==True:
    dfDMEDs={}
    for idx, row in dfcasos.iterrows():
        prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5,"TCSCvar":0.02,"SVCvar":0.02,"UPFCt_sh":0.02,"UPFCV_sh":0.02,"UPFCt_se":0.02,"UPFCV_se":0.05}
        dfDMED=create_DMED(sys,prec,graph,ram,ramUPFC,dfDMEDfp=dDMEDfps[idx])
        dfDMEDFACTs=create_DMED_FACTS(sys,prec,graph,ram,ramUPFC,dfDMEDfp=dDMEDfps[idx])
        dfDMEDsr=pd.concat([dfDMED.copy(),dfDMEDFACTs.copy()])
        dfDMEDs[idx]=dfDMEDsr.copy()
else:
    dfDMEDs={}
    for idx, row in dfcasos.iterrows():
        prec={"SCADAPF":0.02,"SCADAPI":0.02,"SCADAV":0.01,"SMP":0.05,"SMV":0.03,"PSEUDO":0.3,"VIRTUAL":1e-5,"TCSCvar":0.01,"SVCvar":0.01,"UPFCt_sh":0.01,"UPFCV_sh":0.01,"UPFCt_se":0.01,"UPFCV_se":0.01}
        dfDMED=create_DMED(sys,prec,graph,ram,ramUPFC,dfDMEDfp=dDMEDfps[idx])
        dfDMEDs[idx]=dfDMED.copy()


#%%
TCSCini=-0.05
Bini=-0.1
V_sh_ini=1.0
t_sh_ini=0
V_se_ini=0.05
t_se_ini=90*np.pi/180
cx="x2"
dfcasos=pd.DataFrame(data={"TCSC":[20,10,-15,-20],"SVC":[-1,-1,1,1],"UPFC_flow":[20,10,-10,-20],"UPFC_V":[-1,-1,1,1],"TCSC_ini":[TCSCini,TCSCini,TCSCini,TCSCini],"SVC_ini":[Bini,Bini,Bini,Bini]})
conv_LMs={}
conv_BCs={}
conv_noBCs={}
nits_LMs={}
nits_BCs={}
nits_noBCs={}
N=100

dState_LM={}
dStateFACTS_LM={}
dState_BC={}
dStateFACTS_BC={}
dState_noBC={}
dStateFACTS_noBC={}
dfITS=pd.DataFrame()
dconv={}
dconv["n"]=[]
dconv["convLM"]=[]
dconv["convGN"]=[]
dconv["convGNbc"]=[]
dconv["nitsLM"]=[]
dconv["nitsGN"]=[]
dconv["nitsGNbc"]=[]
dconv["caso"]=[]


for idx, row in dfcasos.iterrows():
    conv_LMs[idx]=[]
    conv_BCs[idx]=[]
    conv_noBCs[idx]=[]
    nits_LMs[idx]=[]
    nits_BCs[idx]=[]
    nits_noBCs[idx]=[]
    dState_LM[idx]=[]
    dStateFACTS_LM[idx]=[]
    dState_BC[idx]=[]
    dStateFACTS_BC[idx]=[]
    dState_noBC[idx]=[]
    dStateFACTS_noBC[idx]=[]


    for n in range(N): 
        dfDMED=insert_res(dfDMEDs[idx],n)
        # dfDMED=dfDMEDs[idx].copy()
        for key ,tcsc in ramTCSC.items():
            true=dStateTCSC_ref[idx][(dStateTCSC_ref[idx]["tipo"]=="x_tcsc")&(dStateTCSC_ref[idx]["de"]==key)]["val"].values[0]
            # tcsc.xtcsc_ini=true*(1+cx)
            tcsc.xtcsc_ini=row["TCSC_ini"]
        for key,svc in busSVC.items():
            svc.Bini=row["SVC_ini"]
            # true=dStateTCSC_ref[idx][(dStateTCSC_ref[idx]["tipo"]=="B_svc")&(dStateTCSC_ref[idx]["de"]==key)]["val"].values[0]
            # svc.Bini=true*(1+cx)
        for key,upfc in ramUPFC.items():
            # svc.Bini=row["SVC_ini"]
            # true_VSH=dStateTCSC_ref[idx][(dStateTCSC_ref[idx]["tipo"]=="UPFC_Vsh")&(dStateTCSC_ref[idx]["de"]==key)]["val"].values[0]
            # true_TSH=dStateTCSC_ref[idx][(dStateTCSC_ref[idx]["tipo"]=="UPFC_tsh")&(dStateTCSC_ref[idx]["de"]==key)]["val"].values[0]
            # true_VSE=dStateTCSC_ref[idx][(dStateTCSC_ref[idx]["tipo"]=="UPFC_Vse")&(dStateTCSC_ref[idx]["de"]==key)]["val"].values[0]
            # true_TSE=dStateTCSC_ref[idx][(dStateTCSC_ref[idx]["tipo"]=="UPFC_tse")&(dStateTCSC_ref[idx]["de"]==key)]["val"].values[0]
            # upfc.Vsh_ini=true_VSH*(1+cx)
            # upfc.tsh_ini=true_TSH*(1+cx)
            # upfc.Vse_ini=true_VSE*(1+cx)
            # upfc.tse_ini=true_TSE*(1+cx)
            upfc.Vsh_ini=V_sh_ini
            upfc.tsh_ini=t_sh_ini
            upfc.Vse_ini=V_se_ini
            upfc.tse_ini=t_se_ini
        try:
            conv_LM,nits_LM,dfITsLM=SS_WLS_FACTS_LM_BC(graph,dfDMED,ind_i,printgrad=0,printres=0,pirntits=1,flatstart=2,tol=1e-5,tol2=1e-4)
        except:
            conv_LM=0
            nits_LM=30
            dfITsLM=pd.DataFrame()
        
        if conv_LM==True:
            dState_LM[idx].append(get_state(graph,n))
            dStateFACTS_LM[idx].append(get_state_FACTS(ramTCSC,busSVC,ramUPFC,n))
        
        conv_BC,nits_BC,dfITsGNbc=SS_WLS_FACTS_withBC(graph,dfDMED,ind_i,flatstart=2,tol=1e-5,tol2=1e-4,printgrad=0,printres=0,pirntits=1)
        
        if conv_BC==True:
            dState_BC[idx].append(get_state(graph,n))
            dStateFACTS_BC[idx].append(get_state_FACTS(ramTCSC,busSVC,ramUPFC,n))

        conv_noBC,nits_noBC,dfITsGN=SS_WLS_FACTS_noBC(graph,dfDMED,ind_i,flatstart=2,tol=1e-5,tol2=1e-4,printgrad=0,printres=0,pirntits=1)
        
        if conv_noBC==True:
            dState_noBC[idx].append(get_state(graph,n))
            dStateFACTS_noBC[idx].append(get_state_FACTS(ramTCSC,busSVC,ramUPFC,n))

        
        conv_LMs[idx].append(conv_LM)
        conv_BCs[idx].append(conv_BC)
        conv_noBCs[idx].append(conv_noBC)
        nits_LMs[idx].append(nits_LM)
        nits_BCs[idx].append(nits_BC)
        nits_noBCs[idx].append(nits_noBC)

        dconv["n"].append(n)
        dconv["convLM"].append(conv_LM)
        dconv["convGN"].append(conv_noBC)
        dconv["convGNbc"].append(conv_BC)
        dconv["nitsLM"].append(nits_LM)
        dconv["nitsGN"].append(nits_noBC)
        dconv["nitsGNbc"].append(nits_BC)
        dconv["caso"].append(idx)

        if not dfITsLM.empty:
            dfITsLM["method"]="LM"
            dfITsLM["caso"]=idx
            dfITsLM["n"]=n
        if not dfITsLM.empty:
            dfITsGNbc["method"]="GNbc"
            dfITsGNbc["caso"]=idx
            dfITsGNbc["n"]=n
        if not dfITsLM.empty:
            dfITsGN["method"]="GN"
            dfITsGN["caso"]=idx
            dfITsGN["n"]=n
        dfITS=pd.concat([dfITS,dfITsLM,dfITsGNbc,dfITsGN])
        







#%%
dfSATES=pd.DataFrame()
dfSATES_FACTS=pd.DataFrame()





for key,item in dState_LM.items():
    if len(item)>0:
        df=pd.concat(item)
        df["method"]="LM"
        df["scenario"]=key
        dfSATES=pd.concat([dfSATES,df])
        
for key,item in dStateFACTS_LM.items():
    if len(item)>0:
        df=pd.concat(item)
        df["method"]="LM"
        df["scenario"]=key
        dfSATES_FACTS=pd.concat([dfSATES_FACTS,df])

for key,item in dState_BC.items():
    if len(item)>0:
        df=pd.concat(item)
        df["method"]="GNbc"
        df["scenario"]=key
        dfSATES=pd.concat([dfSATES,df])

for key,item in dStateFACTS_BC.items():
    if len(item)>0:
        df=pd.concat(item)
        df["method"]="GNbc"
        df["scenario"]=key
        dfSATES_FACTS=pd.concat([dfSATES_FACTS,df])


for key,item in dState_noBC.items():
    if len(item)>0:
        df=pd.concat(item)
        df["method"]="GN"
        df["scenario"]=key
        dfSATES=pd.concat([dfSATES,df])


for key,item in dStateFACTS_noBC.items():
    if len(item)>0:
        df=pd.concat(item)
        df["method"]="GN"
        df["scenario"]=key
        dfSATES_FACTS=pd.concat([dfSATES_FACTS,df])



#%%
dfconv=pd.DataFrame(data=dconv)
dfconv.to_csv("ResultadosPSCC/resultados_conv_"+sys+str(cx)+Meas+".csv")
#%%
dfSATES.to_csv("ResultadosPSCC/state_"+sys+str(cx)+Meas+".csv")
dfSATES_FACTS.to_csv("ResultadosPSCC/state_FACTS_"+sys+str(cx)+Meas+".csv")


#%%