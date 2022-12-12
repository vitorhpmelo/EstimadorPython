#%%
import pandas as pd
import matplotlib.pyplot as plt

sistemas=["14","30","118"]



dNormal_5={}
dQR_5={}
dNormal_6={}
dQR_6={}
dNormal_7={}
dQR_7={}
dNormal_8={}
dQR_8={}
dLagran={}
for sys in sistemas:
    dNormal_5[sys]=pd.read_csv("convIEEE"+sys+"Normal5.csv",header=None)
    dNormal_5[sys].columns=["it","normgrad","maxdx"]
    dNormal_6[sys]=pd.read_csv("convIEEE"+sys+"Normal6.csv",header=None)
    dNormal_6[sys].columns=["it","normgrad","maxdx"]
    dNormal_7[sys]=pd.read_csv("convIEEE"+sys+"Normal7.csv",header=None)
    dNormal_7[sys].columns=["it","normgrad","maxdx"]
    dNormal_8[sys]=pd.read_csv("convIEEE"+sys+"Normal8.csv",header=None)
    dNormal_8[sys].columns=["it","normgrad","maxdx"]
    dQR_5[sys]=pd.read_csv("convIEEE"+sys+"QR5.csv",header=None)
    dQR_5[sys].columns=["it","normgrad","maxdx"]
    dQR_6[sys]=pd.read_csv("convIEEE"+sys+"QR6.csv",header=None)
    dQR_6[sys].columns=["it","normgrad","maxdx"]
    dQR_7[sys]=pd.read_csv("convIEEE"+sys+"QR7.csv",header=None)
    dQR_7[sys].columns=["it","normgrad","maxdx"]
    dQR_8[sys]=pd.read_csv("convIEEE"+sys+"QR8.csv",header=None)
    dQR_8[sys].columns=["it","normgrad","maxdx"]
    dLagran[sys]=pd.read_csv("convIEEE"+sys+"lagran.csv",header=None)
    dLagran[sys].columns=["it","normgrad","maxdx"]


fig, ax= plt.subplots(nrows=3,ncols=3,figsize=(14,8))
#%%
i=0
for sys in sistemas:
    ax[0][i].set_title("Equação normal IEEE "+sys,fontsize="x-large")
    normaliz=dNormal_5[sys].normgrad[0]
    ax[0][i].semilogy(dNormal_5[sys].it,dNormal_5[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-5}$",ls="--",marker="d",color="k") 
    normaliz=dNormal_6[sys].normgrad[0]
    ax[0][i].semilogy(dNormal_6[sys].it,dNormal_6[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-6}$",ls="--",marker="o",color="g") 
    normaliz=dNormal_7[sys].normgrad[0]
    ax[0][i].semilogy(dNormal_7[sys].it,dNormal_7[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-7}$",ls="--",marker="x",color="m")    
    normaliz=dNormal_8[sys].normgrad[0]
    ax[0][i].semilogy(dNormal_8[sys].it,dNormal_8[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-8}$",ls="--",marker="*",color="b") 
    ax[0][i].hlines(y=1e-9,xmin=0,xmax=8,colors="r",label="crit de conv")   
    ax[1][i].set_title("Equação QR IEEE "+sys,fontsize="x-large")
    normaliz=dQR_5[sys].normgrad[0]
    ax[1][i].semilogy(dQR_5[sys].it,dQR_5[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-5}$",ls="--",marker="d",color="k") 
    normaliz=dQR_6[sys].normgrad[0]
    ax[1][i].semilogy(dQR_6[sys].it,dQR_6[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-6}$",ls="--",marker="o",color="g") 
    normaliz=dQR_7[sys].normgrad[0]
    ax[1][i].semilogy(dQR_7[sys].it,dQR_7[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-7}$",ls="-.",marker="x",color="m")    
    normaliz=dQR_8[sys].normgrad[0]
    ax[1][i].semilogy(dQR_8[sys].it,dQR_8[sys].normgrad/normaliz,label=r"$\sigma_{virtual} = 10^{-8}$",ls=":",marker="*",color="b")    
    ax[1][i].hlines(y=1e-9,xmin=0,xmax=max(dQR_8[sys].it),colors="r",label="crit de conv")   
    ax[2][i].set_title("Lagrangeano IEEE "+sys,fontsize="x-large")
    normaliz=dLagran[sys].normgrad[0]
    ax[2][i].semilogy(dLagran[sys].it,dLagran[sys].normgrad/normaliz,marker="d",color="k") 
    ax[2][i].hlines(y=1e-6,xmin=0,xmax=max(dLagran[sys].it),colors="r",label="crit de conv")   

    
    ax[0][i].set_ylim(top=1e3)
   
    ax[0][i].set_xticks(dNormal_8[sys].it)
    ax[1][i].set_xticks(dQR_8[sys].it)
    ax[2][i].set_xticks(dLagran[sys].it)
    ax[0][i].set_xlim([0,6])

    ax[0][i].legend()
    ax[1][i].legend()
    ax[2][i].legend()

    ax[0][i].set_xlabel("Iteração",fontsize="x-large")
    ax[1][i].set_xlabel("Iteração",fontsize="x-large")
    ax[2][i].set_xlabel("Iteração",fontsize="x-large")
    s=r"$ \frac{\vert \vert\nabla f(x_k) \vert \vert}{\vert \vert\nabla f(x_0)\vert \vert} $"
    ax[0][i].set_ylabel(s,rotation=0,labelpad=18,fontsize="x-large")
    ax[1][i].set_ylabel(s,rotation=0,labelpad=18,fontsize="x-large")
    ax[2][i].set_ylabel(s,rotation=0,labelpad=18,fontsize="x-large")

    ax[0][i].grid()
    ax[1][i].grid()
    ax[2][i].grid()
    i=i+1
# %%
plt.tight_layout()
fig.savefig("conv.pdf")
fig.savefig("conv.png")
plt.show()
# %%
