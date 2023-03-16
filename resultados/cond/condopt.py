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
k=0.62
fig, ax = plt.subplots(ncols=len(sys),figsize=(11*k,5*k))

xrotulos=["1e-5","1e-6","1e-7","1e-8","1e-9"]
w=0.4
i=0
# fig.suptitle("Condicionamento númerico X Desvio padrão das medidas virtuais",size="x-large")
for system in sys: 
    ax[i].set_title("IEEE "+system)
    
    ax[i].set_xlabel("$k$",fontsize="large")
    ax[i].set_ylabel(r"$\sigma_{virtuais}$",fontsize="large")    
    ax[i].set_yticks(list(range(5)))
    ax[i].set_yticklabels(xrotulos)
    ax[i].barh(list(range(5)),conds[system][5:10],height=w,label="EnQR",color="b")
    ax[i].barh(np.array(list(range(5)))+(w*1.05),conds[system][0:5],height=w,label="WLS-T",color="r")
    ax[i].set_xscale('log')
    ax[i].grid()
    ax[i].legend()
    i=i+1

fig.tight_layout()
fig.savefig("condbar.pdf")
plt.close()