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


dlvls={"5":"+5%","10":"+10%","15":"+15%","20":"+20%","m5":"-5%","m10":"-10%","m15":"-15%","m20":"-20%"}

#%%

teste="Teste3"
sys="IEEE14"

dfcondsAT3_14={}
dfcondsBT3_14={}
dfconvAT3_14={}
dfconvBT3_14={}
for lvl in lvls:
    dfcondsAT3_14[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT3_14[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"B.csv",header=None,skiprows=1)
    dfconvAT3_14[lvl]=pd.read_csv(teste+"/"+sys+"/conv_A_"+lvl+".csv")
    dfconvBT3_14[lvl]=pd.read_csv(teste+"/"+sys+"/conv_B_"+lvl+".csv")

#%%


teste="Teste1"
sys="IEEE14"
dfcondsAT1_14={}
dfcondsBT1_14={}
dfconvAT1_14={}
dfconvBT1_14={}

for lvl in lvls:
    dfcondsAT1_14[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT1_14[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"B.csv",header=None,skiprows=1)
    dfconvAT1_14[lvl]=pd.read_csv(teste+"/"+sys+"/conv_a_"+lvl+".csv")
    dfconvBT1_14[lvl]=pd.read_csv(teste+"/"+sys+"/conv_b_"+lvl+".csv")

#%%




#%%
dfcondsAT3_118={}
dfcondsBT3_118={}
dfconvAT3_118={}
dfconvBT3_118={}

for lvl in lvls:
    dfcondsAT3_118[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT3_118[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"B.csv",header=None,skiprows=1)
    dfconvAT3_118[lvl]=pd.read_csv(teste+"/"+sys+"/conv_a_"+lvl+".csv")
    dfconvBT3_118[lvl]=pd.read_csv(teste+"/"+sys+"/conv_b_"+lvl+".csv")
#%%


teste="Teste1"
sys="IEEE118"
dfcondsAT1_118={}
dfcondsBT1_118={}
dfconvAT1_118={}
dfconvBT1_118={}
for lvl in lvls:
    dfcondsAT1_118[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"A.csv",header=None,skiprows=1)
    dfcondsBT1_118[lvl]=pd.read_csv(teste+"/"+sys+"/cond"+lvl+"B.csv",header=None,skiprows=1)
    dfconvAT1_118[lvl]=pd.read_csv(teste+"/"+sys+"/conv_a_"+lvl+".csv")
    dfconvBT1_118[lvl]=pd.read_csv(teste+"/"+sys+"/conv_b_"+lvl+".csv")

#%%




k=0.7
fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(11*k, 7*k))


#%%

ax[0][1].set_title("Condition number")
ax[0][1].semilogy(range(len(dfcondsAT1_14["5"])),dfcondsAT1_14["5"][0],marker="d",label="FS +5%",color=red)
ax[0][1].semilogy(range(len(dfcondsAT3_14["5"])),dfcondsAT3_14["5"][0],marker="d",label="DC+5%",color=blue)
ax[0][1].semilogy(range(len(dfcondsAT1_14["15"])),dfcondsAT1_14["15"][0],marker="*",label="FS +15%",color=green)
ax[0][1].semilogy(range(len(dfcondsAT3_14["15"])),dfcondsAT3_14["15"][0],marker="*",label="DC+15%",color=orange)
ax[0][1].set_xticklabels(list(range(1,1+len(range(len(dfcondsAT1_14["5"]))))))
ax[0][1].set_xlim(xmin=1,xmax=7)
ax[0][1].set_xlabel("iteration")
ax[0][1].set_ylabel("k")
ax[0][1].legend(fontsize="small")
ax[0][1].grid()



ax[0][0].set_title("Convergence criteria")
ax[0][0].semilogy(range(len(dfconvAT1_14["5"])),dfconvAT1_14["5"]["dz"],marker="d",label="FS+5%",color=red)
ax[0][0].semilogy(range(len(dfconvAT3_14["5"])),dfconvAT3_14["5"]["dz"],marker="d",label="DC+5%",color=blue)
ax[0][0].semilogy(range(len(dfconvAT1_14["15"])),dfconvAT1_14["15"]["dz"],marker="*",label="DC+15%",color=green)
ax[0][0].semilogy(range(len(dfconvAT3_14["15"])),dfconvAT3_14["15"]["dz"],marker="*",label="FS+15%",color=orange)
ax[0][0].set_xticklabels(list(range(1,1+len(range(len(dfconvAT1_14["5"]))))))
ax[0][0].set_xlim(xmin=1,xmax=7)
ax[0][0].set_xlabel("iteration")
s=r"$\frac{\vert \vert \nabla J(x_i) \vert \vert }{\vert\vert \nabla J(x_0) \vert \vert }$"
ax[0][0].set_ylabel(s,fontsize=10)
ax[0][0].legend(fontsize="small")
ax[0][0].grid()



ax[1][1].set_title("Condition number")
ax[1][1].semilogy(range(len(dfcondsAT1_118["5"])),dfcondsAT1_118["5"][0],marker="d",label="FS +5%",color=red)
ax[1][1].semilogy(range(len(dfcondsAT3_118["5"])),dfcondsAT3_118["5"][0],marker="d",label="DC +5%",color=blue)
ax[1][1].semilogy(range(len(dfcondsAT1_118["15"])),dfcondsAT1_118["15"][0],marker="*",label="FS+15%",color=green)
ax[1][1].semilogy(range(len(dfcondsAT3_118["15"])),dfcondsAT3_118["15"][0],marker="*",label="DC +15%",color=orange)
ax[1][1].set_xticklabels(list(range(1,1+len(range(len(dfcondsAT1_118["5"]))))))
ax[1][1].set_xticks(list(range(1,1+len(range(len(dfcondsAT1_118["5"]))))))
ax[1][1].set_xlim(xmin=1,xmax=14)
ax[1][1].set_xlabel("iteration")
ax[1][1].set_ylabel("k")
ax[1][1].legend(fontsize="small")

ax[1][1].grid()



ax[1][0].set_title("Convergence criteria")
ax[1][0].semilogy(range(len(dfconvAT1_118["5"])),dfconvAT1_118["5"]["dz"],marker="d",label="FS +5%",color=red)
ax[1][0].semilogy(range(len(dfconvAT3_118["5"])),dfconvAT3_118["5"]["dz"],marker="d",label="DC +5%",color=blue)
ax[1][0].semilogy(range(len(dfconvAT1_118["15"])),dfconvAT1_118["15"]["dz"],marker="*",label="FS +15%",color=green)
ax[1][0].semilogy(range(len(dfconvAT3_118["15"])),dfconvAT3_118["15"]["dz"],marker="*",label="DC +15%",color=orange)
ax[1][0].set_xticks(list(range(1,1+len(range(len(dfconvAT1_118["5"]))))))
ax[1][0].set_xticklabels(list(range(1,1+len(range(len(dfconvAT1_118["5"]))))))
ax[1][0].set_xlim(xmin=1,xmax=14)
ax[1][0].set_xlabel("iteration")
ax[1][0].set_ylabel(s,fontsize=10)
ax[1][0].legend(fontsize="small")
ax[1][0].grid()

# %%
# plt.show()
fig.tight_layout()
plt.savefig("cond_it.pdf")
# %%
