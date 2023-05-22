
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


red="#da2c38"
blue="#3c1642"
green="#245501"
orange="#ff9914"
pink="#e8ffb7"

clrs=[blue,red,green,orange,pink]



lvls=["20","15","10","5","m5","m10","m15","m20"]


dlvls={"5":"5%","10":"10%","15":"15%","20":"20%","m5":"-5%","m10":"-10%","m15":"-15%","m20":"-20%"}


#%%



teste="Teste1"
sys="IEEE14"
dfcondsAT1_14={}
dfcondsBT1_14={}

for lvl in lvls:
    dfcondsAT1_14[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT1_14[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"B.csv",header=None,skiprows=1)

#%%
teste="Teste3"
sys="IEEE14"

dfcondsAT3_14={}
dfcondsBT3_14={}

for lvl in lvls:
    dfcondsAT3_14[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT3_14[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"B.csv",header=None,skiprows=1)
#%%


Yt1A_14={}
Yt1B_14={}
Yt3A_14={}
Yt3B_14={}
for lvl in lvls:
    Yt1A_14[dlvls[lvl]]=dfcondsAT1_14[lvl][0][0]
    Yt1B_14[dlvls[lvl]]=dfcondsBT1_14[lvl][0][0]
    Yt3A_14[dlvls[lvl]]=dfcondsAT3_14[lvl][0][0]
    Yt3B_14[dlvls[lvl]]=dfcondsBT3_14[lvl][0][0]

#%%
# teste="Teste3"
dfcondsAT3_118={}
dfcondsBT3_118={}

for lvl in lvls:
    dfcondsAT3_118[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT3_118[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"B.csv",header=None,skiprows=1)
#%%


teste="Teste1"
sys="IEEE118"
dfcondsAT1_118={}
dfcondsBT1_118={}

for lvl in lvls:
    dfcondsAT1_118[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT1_118[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"B.csv",header=None,skiprows=1)

#%%


Yt1A_118={}
Yt1B_118={}
Yt3A_118={}
Yt3B_118={}
for lvl in lvls:
    Yt1A_118[dlvls[lvl]]=dfcondsAT1_118[lvl][0][0]
    Yt1B_118[dlvls[lvl]]=dfcondsBT1_118[lvl][0][0]
    Yt3A_118[dlvls[lvl]]=dfcondsAT3_118[lvl][0][0]
    Yt3B_118[dlvls[lvl]]=dfcondsBT3_118[lvl][0][0]


#%%

x14=np.array(list(range(1,len(Yt1A_14)+1)))
x118=np.array(list(range(1,len(Yt1A_118)+1)))
wd=0.2


#%%
k=0.65
fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(12*k, 4.5*k))

ax[0].set_title("IEEE 14 bus")
ax[0].set_yscale('log')
ax[0].set_xlabel("Power Flow Variation")
ax[0].set_ylabel("k")

ax[0].bar(x14-wd*3/2,Yt1A_14.values(),width=wd,align="center",label="flat-start A",color=clrs[0])
ax[0].bar(x14-wd/2,Yt1B_14.values(),width=wd,align="center",label="flat-start B",color=clrs[1])
ax[0].bar(x14+wd/2,Yt3A_14.values(),width=wd,align="center",label="DC start A",color=clrs[2])
ax[0].bar(x14+wd*3/2,Yt3B_14.values(),width=wd,align="center",label="DC start B",color=clrs[3])
ax[0].set_xticks(x14)
ax[0].set_xticklabels(list(Yt1A_14.keys()),rotation=0)
ax[0].legend()
ax[0].grid()

ax[1].set_title("IEEE 118 bus")
ax[1].set_yscale('log')
ax[1].set_xlabel("Power Flow Variation")
ax[1].set_ylabel("k")
ax[1].bar(x118-wd*3/2,Yt1A_118.values(),width=wd,align="center",label="flat-start A",color=clrs[0])
ax[1].bar(x118-wd/2,Yt1B_118.values(),width=wd,align="center",label="flat-start B",color=clrs[1])
ax[1].bar(x118+wd/2,Yt3A_118.values(),width=wd,align="center",label="DC start A",color=clrs[2])
ax[1].bar(x118+wd*3/2,Yt3B_118.values(),width=wd,align="center",label="DC start B",color=clrs[3])
ax[1].set_xticks(x118)
ax[1].set_xticklabels(list(Yt1A_118.keys()),rotation=0)
ax[1].legend()
ax[1].grid()

plt.tight_layout()

fig.savefig("cond_teste1.pdf")

# %%
