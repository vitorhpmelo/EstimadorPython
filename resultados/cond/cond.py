#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Paleta de cores do Artigo SBSE



color1="#391DC2"
color2="#BF1E17"
color3="#188C6C"
color4="#FBC519"



#%%



sys=["30","118"]
conds={}

for system in sys:
    conds[system]=np.loadtxt("condIEEE"+system+".csv",dtype=float,delimiter=",")
#%%
k=0.62
fig, ax = plt.subplots(ncols=len(sys),figsize=(10*k,4*k))

xrotulos=["1e-5","1e-6","1e-7","1e-8","1e-9"]

i=0
# fig.suptitle("Condicionamento númerico X Desvio padrão das medidas virtuais",size="x-large")
for system in sys: 
    ax[i].set_title("IEEE "+system)
    ax[i].semilogy(list(range(5)),conds[system][0:5],marker="d",label="WLS-T",color=color1)
    ax[i].set_ylabel("$k$",fontsize="large")
    ax[i].set_xlabel(r"$\sigma_{virtuais}$",fontsize="large")    
    ax[i].set_xticks(list(range(5)))
    ax[i].set_xticklabels(xrotulos)
    ax[i].semilogy(list(range(5)),conds[system][5:10],marker="o",label="EnQR",color=color2)
    ax[i].grid()
    ax[i].legend()
    i=i+1

fig.tight_layout()
fig.savefig("cond.pdf")
plt.close()
# %%
