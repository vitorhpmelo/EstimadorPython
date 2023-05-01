#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


red="#da2c38"
blue="#3c1642"
green="#245501"
orange="#ff9914"
pink="#e8ffb7"

clrs=[blue,red,green,orange,pink]

lvls=["20","15","10","5"]


dlvls={"5":"+5%","10":"+10%","15":"+15%","20":"+20%","m5":"-5%","m10":"-10%","m15":"-15%","m20":"-20%"}

dsys={"IEEE 14":"IEEE14_tcsc","IEEE 118":"IEEE118_tcsc_2"}


#%%
path="EMA/"

df={}
for key,item in dsys.items():
    df[key]=pd.DataFrame()
    for lvl in lvls:
        d=pd.read_csv(path+"res"+item+"str_5incr_"+str(lvl)+".csv")
        df[key]=pd.concat([df[key],d])
# %%

k=0.60
x=8
y=6.5

fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(x*k, y*k))


i=0
wd=0.2

for key,sys in dsys.items():
    x=np.array(list(range(1,len(df[key][(df[key]["Solver"]=="A")])+1)))
    axs[i].set_ylabel("MAE")
    axs[i].set_xlabel("Power Flow Increase %")
    axs[i].set_xticks(range(1,len(df[key][(df[key]["Solver"]=="A")]["prec"])+1))
    axs[i].set_xticklabels(df[key][(df[key]["Solver"]=="A")]["prec"].tolist())
    axs[i].set_title(key)
    axs[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axs[i].bar(x-wd,df[key][(df[key]["Solver"]=="A")]["EMA_V"],width=wd,label="MAE V",color=red)
    axs[i].bar(x,df[key][(df[key]["Solver"]=="A")]["EMA_T"],width=wd,label="MAE $\\theta$",color=blue)
    axs[i].bar(x+wd,df[key][(df[key]["Solver"]=="A")]["EMA_X"],width=wd,label="MAE X",color=green)
    axs[i].grid()
    axs[i].legend()
    i=i+1
    
# %%
fig.tight_layout()
plt.savefig("EMA.pdf")
# %%
