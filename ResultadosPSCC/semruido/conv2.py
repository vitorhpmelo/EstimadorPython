
#%%
# -*- coding: utf-8 -*-



import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np


dsys={"IEEE118":"118 bus"}
colors=["#D6870D","#099358","#D60D0D","#291C95"]
dmeas={"":"FACTS not measured","MEAS":"FACTS Measured"}
dfDATA={}
sys="IEEE118"
#%%
for key, name in dmeas.items():
    dfx1=pd.read_csv("Its_"+sys+"_rakp2009"+key+"x1.csv")
    dfx2=pd.read_csv("Its_"+sys+"_rakp2009"+key+"x2.csv")
    dfx1["ini"]=1
    dfx2["ini"]=2
    dfDATA[key]=pd.concat([dfx1,dfx2])


# %%
k=0.55
fig, ax =plt.subplots(ncols=2,nrows=1,figsize=(k*15,k*6))
6
caso=1
i=0
axprime=ax.copy()
d={0:"a) ",1:"b) "}
fstitle=18
fslabel=15
fslegend=12

for key, name in dmeas.items():
    xmax=-1
    ax[i].set_title(d[i]+name,fontsize=fstitle)
    mask1=(dfDATA[key]["caso"]==caso)&(dfDATA[key]["ini"]==1)
    mask2=(dfDATA[key]["method"]=="LM")

    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color=colors[1],label="LM",marker="o")
    xmax=len(dfDATA[key][mask1&mask2])
    
    mask1=(dfDATA[key]["caso"]==caso)&(dfDATA[key]["ini"]==1)
    mask2=(dfDATA[key]["method"]=="GN")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color=colors[2],label="GN",marker="o")
    xmax=max(xmax,len(dfDATA[key][mask1&mask2]))
    mask1=(dfDATA[key]["caso"]==caso)&(dfDATA[key]["ini"]==1)
    mask2=(dfDATA[key]["method"]=="GNbc")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color=colors[3],label="GNbc",marker="o")
    xmax=max(xmax,len(dfDATA[key][mask1&mask2]))
    mask1=(dfDATA[key]["caso"]==caso)&(dfDATA[key]["ini"]==2)
    mask2=(dfDATA[key]["method"]=="LM")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color=colors[1],label="LM",ls="--",marker="x")
    xmax=max(xmax,len(dfDATA[key][mask1&mask2]))
    mask1=(dfDATA[key]["caso"]==caso)&(dfDATA[key]["ini"]==2)
    mask2=(dfDATA[key]["method"]=="GN")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color=colors[2],label="GN",ls="--",marker="x")
    xmax=max(xmax,len(dfDATA[key][mask1&mask2]))
    mask1=(dfDATA[key]["caso"]==caso)&(dfDATA[key]["ini"]==2)
    mask2=(dfDATA[key]["method"]=="GNbc")
    ax[i].semilogy(range(len(dfDATA[key][mask1&mask2])),dfDATA[key][(mask1)&(mask2)]["dx"],color=colors[3],label="GNbc",ls="--",marker="x")
    xmax=max(xmax,len(dfDATA[key][mask1&mask2]))
    ax[i].set_xlabel("Iteration",fontsize=fslabel)
    ax[i].set_ylabel(r"$\vert \vert \Delta x \vert \vert$",fontsize=fslabel)
   
    ax[i].set_xlim(xmin=0,xmax=xmax)
    ax[i].set_ylim(ymin=1e-5,ymax=1e2)
    ax[i].set_xticks(range(0,xmax,2))

    ax[i].grid()
    LM_patch=mpatches.Patch(color=colors[1],label="LM")
    
    GN_patch=mpatches.Patch(color=colors[2],label="GN")
    GNbc_patch=mpatches.Patch(color=colors[3],label="GNbc")
    x1line=mlines.Line2D([],[],color="k",marker="o",label=r"$x^1_0$")
    x2line=mlines.Line2D([],[],color="k",marker="x",ls="--",label=r"$x^2_0$")
    ax[i].legend(handles=[LM_patch,GNbc_patch,GN_patch],loc=3,title=r"$\bf{Method}$",fancybox=False,fontsize=fslegend)
    
    axprime[i]=ax[i].twinx()
    axprime[i].legend(handles=[x1line,x2line],loc=1,title=r"$\bf{Ini}$",fancybox=False,fontsize=fslegend)
    ax[i].get_shared_y_axes().join(ax[i], axprime[i])
    axprime[i].axes.get_yaxis().set_visible(False)
    i=i+1
# %%




plt.tight_layout()
plt.savefig("conv_semRUIDO_"+sys+".pdf")
#%%