#%%
import pandas as pd
import numpy as np

sys="teste118"
dfbran=pd.read_csv(sys+"bran.csv",header=None)
dfbran.columns=["fbus","tbus","r","x","b","rateA","rateB","rateC","ratio","angle","status","angmin","angmax"]
dfbus=pd.read_csv(sys+"bus.csv",header=None)
dfbus.columns=["bus_i","type","Pd","Qd","Gs","Bs","area","Vm","Va","baseKV","zone","Vmax","Vmin"]
dfgen=pd.read_csv(sys+"gen.csv",header=None)
dfgen.columns=["bus","Pg","Qg","Qmax","Qmin","Vg","mBase","status","Pmax","Pmin","Pc1","Pc2","Qc1min","Qc1max","Qc2min","Qc2max","ramp_agc","ramp_10","ramp_30","ramp_q","apf7"]
# %%

dfDBAR=pd.DataFrame()
dfDBAR["id"]=dfbus["bus_i"]
dfDBAR["type"]=3-dfbus["type"]
dfDBAR["V"]=dfbus["Vm"]
dfDBAR["teta"]=dfbus["Va"]
dfDBAR["Pg"]=np.zeros(len(dfDBAR))
dfDBAR["Qg"]=np.zeros(len(dfDBAR))
dfDBAR["Pd"]=dfbus["Pd"]
dfDBAR["Qd"]=dfbus["Qd"]
dfDBAR["Bs"]=dfbus["Bs"]
#%%
for idx, row in dfgen.iterrows():
    dfDBAR.loc[dfDBAR["id"]==row["bus"],"Pg"]=row["Pg"]
    dfDBAR.loc[dfDBAR["id"]==row["bus"],"Qg"]=row["Qg"]
# %%


dfDBAR.to_csv("DBAR.csv",header=None,index=None)

#%%
dfBRAN=pd.DataFrame()
dfBRAN["id"]=list(range(1,len(dfbran)+1))
dfBRAN["type"]=np.ones(len(dfbran),dtype=int)
trafosidx=dfbran[(dfbran["ratio"]!=0)&(dfbran["ratio"]!=1)].index
dfBRAN.loc[trafosidx,"type"]=2
dfBRAN["from"]=dfbran["fbus"]
dfBRAN["to"]=dfbran["tbus"]
dfBRAN["r"]=dfbran["r"]
dfBRAN["x"]=dfbran["x"]
dfBRAN["b"]=dfbran["b"]
dfBRAN["tap"]=dfbran["ratio"]
dfBRAN.loc[dfBRAN["tap"]==0,"tap"]=1
dfBRAN.to_csv("DBRAN.csv",header=None,index=None)
# %%
