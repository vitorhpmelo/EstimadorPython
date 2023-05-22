
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


IEEE14V_A={}
IEEE14T_A={}
IEEE14X_A={}
IEEE14V_B={}
IEEE14T_B={}
IEEE14X_B={}

sys="IEEE 14"

for lvl,name in dlvls.items():
    IEEE14V_A[lvl]=dMAE[sys][lvl]["EMA_V"][0]
    IEEE14T_A[lvl]=dMAE[sys][lvl]["EMA_T"][0]
    IEEE14X_A[lvl]=dMAE[sys][lvl]["EMA_X"][0]
    IEEE14V_B[lvl]=dMAE[sys][lvl]["EMA_V"][1]
    IEEE14T_B[lvl]=dMAE[sys][lvl]["EMA_T"][1]
    IEEE14X_B[lvl]=dMAE[sys][lvl]["EMA_X"][1]


# %%

IEEE118V_A={}
IEEE118T_A={}
IEEE118X_A={}
IEEE118V_B={}
IEEE118T_B={}
IEEE118X_B={}

sys="IEEE 118"

for lvl,name in dlvls.items():
    IEEE118V_A[lvl]=dMAE[sys][lvl]["EMA_V"][0]
    IEEE118T_A[lvl]=dMAE[sys][lvl]["EMA_T"][0]
    IEEE118X_A[lvl]=dMAE[sys][lvl]["EMA_X"][0]
    IEEE118V_B[lvl]=dMAE[sys][lvl]["EMA_V"][1]
    IEEE118T_B[lvl]=dMAE[sys][lvl]["EMA_T"][1]
    IEEE118X_B[lvl]=dMAE[sys][lvl]["EMA_X"][1]
# %%

k=0.55
fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(12*k, 8*k))



x=np.array(range(1,1+len(IEEE14T_A.values())))
wd=0.10

ax[0][0].set_title("IEEE 14 MAE SE A")
ax[0][0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[0][0].bar(x-1.1*wd,IEEE14V_A.values(),width=wd,align="center",label=r"$V$",color=red)
ax[0][0].bar(x,IEEE14T_A.values(),width=wd,align="center",label=r"$\theta$",color=orange)
ax[0][0].bar(x+1.1*wd,IEEE14X_A.values(),width=wd,align="center",label=r"$X_{TCSC}$",color=blue)
ax[0][0].set_xlabel("Power Flow Increase",fontsize="large")
ax[0][0].set_ylabel("MAE")
ax[0][0].set_xticks(x)
ax[0][0].set_xticklabels(dlvls.values())
ax[0][0].grid()
ax[0][1].set_title("IEEE 14 MAE SE B")
ax[0][1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[0][1].bar(x-1.1*wd,IEEE14V_B.values(),width=wd,align="center",label=r"$V$",color=red)
ax[0][1].bar(x,IEEE14T_B.values(),width=wd,align="center",label=r"$\theta$",color=orange)
ax[0][1].bar(x+1.1*wd,IEEE14X_B.values(),width=wd,align="center",label=r"$X_{TCSC}$",color=blue)
ax[0][1].set_xlabel("Power Flow Increase",fontsize="large")
ax[0][1].set_ylabel("MAE")
ax[0][1].set_xticks(x)
ax[0][1].set_xticklabels(dlvls.values())
ax[0][1].grid()
ax[1][0].set_title("IEEE 118 MAE SE A")
ax[1][0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[1][0].bar(x-1.1*wd,IEEE118V_A.values(),width=wd,align="center",label=r"$V$",color=red)
ax[1][0].bar(x,IEEE118T_A.values(),width=wd,align="center",label=r"$\theta$",color=orange)
ax[1][0].bar(x+1.1*wd,IEEE118X_A.values(),width=wd,align="center",label=r"$X_{TCSC}$",color=blue)
ax[1][0].set_xlabel("Power Flow Increase",fontsize="large")
ax[1][0].set_ylabel("MAE")
ax[1][0].set_xticks(x)
ax[1][0].set_xticklabels(dlvls.values())
ax[1][0].grid()
ax[1][1].set_title("IEEE 118 MAE SE B")
ax[1][1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[1][1].bar(x-1.1*wd,IEEE118V_B.values(),width=wd,align="center",label=r"$V$",color=red)
ax[1][1].bar(x,IEEE118T_B.values(),width=wd,align="center",label=r"$\theta$",color=orange)
ax[1][1].bar(x+1.1*wd,IEEE118X_B.values(),width=wd,align="center",label=r"$X_{TCSC}$",color=blue)
ax[1][1].set_xlabel("Power Flow Increase",fontsize="large")
ax[1][1].set_ylabel("MAE")
ax[1][1].set_xticks(x)
ax[1][1].set_xticklabels(dlvls.values())

ax[1][1].grid()
ax[0][0].legend()
ax[0][1].legend()
ax[1][0].legend()
ax[1][1].legend()

# %%

plt.tight_layout()
plt.savefig("MAE2.pdf")

#%%