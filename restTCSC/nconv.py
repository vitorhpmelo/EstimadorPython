
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


dlvls={"5":"+5%","10":"+10%","15":"+15%","20":"+20%"}

dsys={"IEEE 14":"IEEE14_tcsc","IEEE 118":"IEEE118_tcsc_2"}



#%%
dMAE={}

for sys,file in dsys.items():
    dMAE[sys]={}
    for lvl,name in dlvls.items():
        dMAE[sys][lvl]=pd.read_csv("res"+file+"str_5"+"incr_"+lvl+".csv")

#%%


IEEE14_A={}
IEEE14_B={}

sys="IEEE 14"

for lvl,name in dlvls.items():
    IEEE14_A[lvl]=dMAE[sys][lvl]["Namos"][0]-1
    IEEE14_B[lvl]=dMAE[sys][lvl]["Namos"][1]-1


# %%

IEEE118_A={}
IEEE118_B={}

sys="IEEE 118"

for lvl,name in dlvls.items():
    IEEE118_A[lvl]=dMAE[sys][lvl]["Namos"][0]-1
    IEEE118_B[lvl]=dMAE[sys][lvl]["Namos"][1]-1

# %%

k=0.55
fig,ax1 = plt.subplots(figsize=(12*k, 4*k))



x=np.array(range(1,5))
wd=0.10

ax1.set_title("Number of Monte Carlo samples")
ax1.bar(x-3*1.1*wd/2,IEEE14_A.values(),width=wd,align="center",label="IEEE 14 A",color=red)
ax1.bar(x-1.1*wd/2,IEEE14_B.values(),width=wd,align="center",label="IEEE 14 B",color=orange)
ax1.bar(x+1.1*wd/2,IEEE118_A.values(),width=wd,align="center",label="IEEE 118 A",color=blue)
ax1.bar(x+3*1.1*wd/2,IEEE118_B.values(),width=wd,align="center",label="IEEE 118 B",color=green)
ax1.set_xlabel("Power Flow Increase",fontsize="large")
ax1.set_ylabel("Number")
ax1.set_xticks(x)
ax1.set_xticklabels(dlvls.values())
ax1.grid()
ax1.legend(ncols=2)
ax1.set_ylim(ymin=80)
# %%

plt.tight_layout()
plt.savefig("nsamples.pdf")

#%%