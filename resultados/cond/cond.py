#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%

sys=["14","30","118"]
conds={}

for system in sys:
    conds[system]=np.loadtxt("condIEEE"+system+".csv",dtype=float,delimiter=",")
#%%
fig, ax = plt.subplots(ncols=3,figsize=(14,4))

xrotulos=["1e-5","1e-6","1e-7","1e-8","1e-9"]

#%%
i=0
fig.suptitle("Condicionamento númerico X Desvio padrão das medidas virtuais")
for system in sys: 
    ax[i].set_title("Sistema IEEE "+system)
    ax[i].semilogy(list(range(5)),conds[system][0:5],marker="d",label="Normal",color="r")
    ax[i].set_ylabel("Número de condicionamento")
    ax[i].set_xlabel("Desvio Padrão ")    
    ax[i].set_xticks(list(range(5)))
    ax[i].set_xticklabels(xrotulos)
    ax[i].semilogy(list(range(5)),conds[system][5:10],marker="o",label="QR",color="b")
    ax[i].grid()
    ax[i].legend()
    i=i+1

plt.tight_layout()
fig.savefig("cond.pdf")
# %%
