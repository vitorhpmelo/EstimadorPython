
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




teste="Teste3"
sys="IEEE14"

dfcondsAT3_14={}


for lvl in lvls:
    dfcondsAT3_14[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"A.csv",header=None,skiprows=1)
#%%


Yt3A_14={}
for lvl in lvls:
    Yt3A_14[dlvls[lvl]]=np.array(dfcondsAT3_14[lvl][0])[-1]

#%%
# teste="Teste3"
dfcondsAT3_118={}

teste="Teste3"
sys="IEEE118"

for lvl in lvls:
    dfcondsAT3_118[lvl]=pd.read_csv(teste+"/"+sys+"/conds_"+lvl+"A.csv",header=None,skiprows=1)
#%%


Yt3A_118={}
Yt3B_118={}
for lvl in lvls:
    Yt3A_118[dlvls[lvl]]=np.array(dfcondsAT3_118[lvl][0])[-1]


#%%

x14=np.array(list(range(1,len(lvls)+1)))
x118=np.array(list(range(1,len(lvls)+1)))
wd=0.2


#%%
k=0.55
fig,ax1 = plt.subplots(figsize=(12*k, 4*k))
ax2 = ax1.twinx()

ax1.set_title("Condition number at the last iteration")
# ax1.set_yscale('log')
ax1.set_xlabel("Power Flow Variation")
ax1.set_ylabel("k IEEE 14", color=clrs[0])

ax1.bar(x14-wd/2,Yt3A_14.values(),width=wd,align="center",label="IEEE 14",color=clrs[0])
ax1.set_xticks(x14)
ax1.set_xticklabels(list(Yt3A_14.keys()),rotation=0)
ax1.legend()
ax1.grid()


ax2.set_ylabel("k IEEE 118", color=clrs[3])
ax2.bar(x14+wd/2,Yt3A_118.values(),width=wd,align="center",label="IEEE 118",color=clrs[3])

ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
plt.tight_layout()

fig.savefig("cond_teste2.pdf")

# %%
