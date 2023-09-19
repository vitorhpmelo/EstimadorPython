#%%
# -*- coding: utf-8 -*-



import pandas as pd
import matplotlib.pyplot as plt


dsys={"IEEE14":"14 bus","IEEE118":"118 bus"}


dfDATA={}

for key, sys in dsys.items():
    dfx1=pd.read_csv(key+"/Its_x1.csv")
    dfx2=pd.read_csv(key+"/Its_x2.csv")
    dfx1["ini"]=1
    dfx2["ini"]=2
    dfDATA[key]=pd.concat([dfx1,dfx2])


# %%

fig, ax =plt.subplots(ncols=2,nrows=1,figsize=(16,16))


i=0
for key, sys in dsys.items():
    mask1=(dfDATA[key]["caso"]==1)&(dfDATA[key]["ini"]==1)
    mask2=(dfDATA[key]["method"]=="LM")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color="b",label="LM")
    mask1=(dfDATA[key]["caso"]==1)&(dfDATA[key]["ini"]==1)
    mask2=(dfDATA[key]["method"]=="GN")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color="r",label="GN")
    mask1=(dfDATA[key]["caso"]==1)&(dfDATA[key]["ini"]==1)
    mask2=(dfDATA[key]["method"]=="GNbc")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color="g",label="GNbc")
    
    i=i+1
# %%
plt.show()