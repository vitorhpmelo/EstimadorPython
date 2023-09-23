#%%
import pandas as pd
import numpy as np
import numpy.linalg as liang
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['pgf.texsystem'] = 'pdflatex'
mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

colors=["#D6870D","#099358","#D60D0D","#291C95"]



dsys={"IEEE14_rakp2009":"14 bus","IEEE118_rakp2009":"118 bus"}




dfDATA=pd.read_csv("Resultado_conv_ruido.csv") 
# %%

k=0.55

fig, ax =plt.subplots(ncols=2,nrows=1,figsize=(14*k,6*k))



wd=0.20
off=wd/8
i=0
AB=["a)","b)"]

fslabels=16
fstitle=16
fslegend=13
dtotalconv={}
for key, sys in dsys.items():
    
    dtotalconv[sys]={}
    
    ylm_x1=dfDATA[(dfDATA["sys"]==key) & (dfDATA["ini"]==1)&(dfDATA["Med"]==0)]["convLM"].array
    dtotalconv[sys]["LM"]=sum(ylm_x1)
    yGN_x1=dfDATA[(dfDATA["sys"]==key) & (dfDATA["ini"]==1)&(dfDATA["Med"]==0)]["convGN"].array
    dtotalconv[sys]["GN"]=sum(yGN_x1)
    yGNbc_x1=dfDATA[(dfDATA["sys"]==key) & (dfDATA["ini"]==1)&(dfDATA["Med"]==0)]["convGNbc"].array
    dtotalconv[sys]["GNbc"]=sum(yGNbc_x1)
    
    x1=np.arange(0,len(ylm_x1))*2+1

    ax[i].set_title(AB[i]+sys,fontsize=fstitle)
    ax[i].bar(x1+(wd+off),yGN_x1,width=wd,color=colors[2],label="GN")
    ax[i].bar(x1,ylm_x1,width=wd,color=colors[1],label="LM")
    ax[i].bar(x1-(wd+off),yGNbc_x1,width=wd,color=colors[3],label="GNbc")
    ax[i].grid(axis="y")
    ax[i].set_ylabel("nº convs",fontsize=fslabels)
    ax[i].set_xlabel("Secenarios/Ini",fontsize=fslabels)
    

    ylm_x2=dfDATA[(dfDATA["sys"]==key) & (dfDATA["ini"]==2)&(dfDATA["Med"]==0)]["convLM"].array
    dtotalconv[sys]["LM"]=dtotalconv[sys]["LM"]+sum(ylm_x2)
    yGN_x2=dfDATA[(dfDATA["sys"]==key) & (dfDATA["ini"]==2)&(dfDATA["Med"]==0)]["convGN"].array
    dtotalconv[sys]["GN"]=dtotalconv[sys]["GN"]+sum(yGN_x2)
    yGNbc_x2=dfDATA[(dfDATA["sys"]==key) & (dfDATA["ini"]==2)&(dfDATA["Med"]==0)]["convGNbc"].array
    dtotalconv[sys]["GNbc"]=dtotalconv[sys]["GNbc"]+sum(yGNbc_x2)
    x2=np.arange(0,len(ylm_x1))*2+2
    ax[i].bar(x2+(wd+off),yGN_x2,width=wd,color=colors[2])
    ax[i].bar(x2,ylm_x2,width=wd,color=colors[1])
    ax[i].bar(x2-(wd+off),yGNbc_x2,width=wd,color=colors[3])

    xticks=np.sort(np.concatenate([x1,x2]))
    ax[i].set_xticks(xticks)
    ax[i].legend(loc="lower left",fontsize=fslegend)
    i=i+1
# %%

lst=["A","B","C","D"]
i=0
xtickslabels=list(xticks.copy())
x1labels=[]
for x in x1:
    s=r'$\begin{bmatrix}'+lst[i] + r'\\ x_0^1 \end{bmatrix}$'
    xtickslabels[x-1]=s
    i=i+1
i=0
for x in x2:
    s=r'$\begin{bmatrix}'+lst[i] + r'\\ x_0^2 \end{bmatrix}$'
    xtickslabels[x-1]=s
    i=i+1

i=0

#%%
for key, sys in dsys.items():
    ax[i].set_xticklabels(xtickslabels,fontsize=fslabels)
    ax[i].tick_params(axis="y", rotation = 90) 

    ax[i].set_yticks([0,25,50,75,100])    
    ax[i].set_yticklabels([0,25,50,75,100],fontsize=fslegend) 
    i=i+1

# %%
plt.tight_layout()
plt.savefig("nconvs_sem_meds.pdf")
# %%
