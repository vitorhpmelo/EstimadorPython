
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
dMAE_5={}


for sys,file in dsys.items():
    dMAE_5[sys]={}
    for lvl,name in dlvls.items():
        dMAE_5[sys][lvl]=pd.read_csv("res"+file+"str_5"+"incr_"+lvl+".csv")

dMAE_4={}
for sys,file in dsys.items():
    dMAE_4[sys]={}
    for lvl,name in dlvls.items():
        dMAE_4[sys][lvl]=pd.read_csv("res"+file+"str_4"+"incr_"+lvl+".csv")

#%%


IEEE14Total_A_5={}
IEEE14Total_B_5={}
IEEE14Total_A_4={}
IEEE14Total_B_4={}

sys="IEEE 14"

for lvl,name in dlvls.items():
    IEEE14Total_A_5[lvl]=dMAE_5[sys][lvl]["EMA_total"][0]
    IEEE14Total_B_5[lvl]=dMAE_5[sys][lvl]["EMA_total"][1]
    IEEE14Total_A_4[lvl]=dMAE_4[sys][lvl]["EMA_total"][0]
    IEEE14Total_B_4[lvl]=dMAE_4[sys][lvl]["EMA_total"][1]



# %%

IEEE118Total_A_5={}
IEEE118Total_B_5={}
IEEE118Total_A_4={}
IEEE118Total_B_4={}

sys="IEEE 118"

for lvl,name in dlvls.items():
    IEEE118Total_A_5[lvl]=dMAE_5[sys][lvl]["EMA_total"][0]
    IEEE118Total_B_5[lvl]=dMAE_5[sys][lvl]["EMA_total"][1]
    IEEE118Total_A_4[lvl]=dMAE_4[sys][lvl]["EMA_total"][0]
    IEEE118Total_B_4[lvl]=dMAE_4[sys][lvl]["EMA_total"][1]

# %%

k=0.55
fig,ax = plt.subplots(ncols=2,nrows=1,figsize=(12*k, 4*k))



x=np.array(range(0,len(IEEE14Total_A_5.values())))
wd=0.10
#%%

x14_A_4=x[~np.isnan(np.array(list(IEEE14Total_A_4.values())))]
x14_B_4=x[~np.isnan(np.array(list(IEEE14Total_B_4.values())))]
#%%
ax[0].set_title("IEEE 14 MAE")
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[0].bar(x-1.1*wd/2,IEEE14Total_A_5.values(),width=wd,align="center",label=r"SE A DC",color=red)
ax[0].bar(x-3*1.1*wd/2,IEEE14Total_B_5.values(),width=wd,align="center",label=r"SE-B DC",color=blue)
y=np.array(list(IEEE14Total_A_4.values()))[x14_A_4]
ax[0].bar(x14_A_4+1.1*wd/2,y,width=wd,align="center",label=r"SE A FS",color=orange)
y=np.array(list(IEEE14Total_B_4.values()))[x14_B_4]
ax[0].bar(x14_B_4+3*1.1*wd/2,y,width=wd,align="center",label=r"SE B FS",color=green)
ax[0].set_xlabel("Power Flow Increase",fontsize="large")
ax[0].set_ylabel("MAE")
ax[0].set_xticks(x)
ax[0].set_xticklabels(dlvls.values())
ax[0].grid()

#%%



x118_A_4=x[~np.isnan(np.array(list(IEEE118Total_A_4.values())))]
x118_B_4=x[~np.isnan(np.array(list(IEEE118Total_B_4.values())))]

#%%
ax[1].set_title("IEEE 118 MAE")
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[1].bar(x-1.1*wd/2,IEEE118Total_A_5.values(),width=wd,align="center",label=r"SE A DC",color=red)
ax[1].bar(x-3*1.1*wd/2,IEEE118Total_B_5.values(),width=wd,align="center",label=r"SE B DC",color=blue)
y=np.array(list(IEEE118Total_A_4.values()))[x118_A_4]
ax[1].bar(x118_A_4+1.1*wd/2,y,width=wd,align="center",label=r"SE A FS",color=orange)
y=np.array(list(IEEE118Total_B_4.values()))[x118_B_4]
ax[1].bar(x118_B_4+3*1.1*wd/2,y,width=wd,align="center",label=r"SE B FS",color=green)
ax[1].set_xlabel("Power Flow Increase",fontsize="large")
ax[1].set_ylabel("MAE")
ax[1].set_xticks(x)
ax[1].set_xticklabels(dlvls.values())

ax[1].grid()
ax[0].legend()

ax[1].legend()


# %%

plt.tight_layout()
plt.savefig("MAE3.pdf")

#%%