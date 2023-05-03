
#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#%%
lstsys=["IEEE14","IEEE118"]

G={}
Lagr={}
Rqr={}
for sys in lstsys:
    G[sys]=np.loadtxt(sys+"/G.csv",dtype=float,delimiter=",")
    Lagr[sys]=np.loadtxt(sys+"/Lagra.csv",dtype=float,delimiter=",")
    Rqr[sys]=np.loadtxt(sys+"/Rqr.csv",dtype=float,delimiter=",")

k=0.60
fig, ax =plt.subplots(ncols=4,figsize=(k*15,k*6))
#%%
nnzG={}
nnzLag={}
nnzRqr={}
nnzGpercentual={}
nnzLagpercentual={}
nnzRqrpercentual={}
for sys in lstsys:
    nnzG[sys]=np.count_nonzero(G[sys])
    nnzLag[sys]=np.count_nonzero(Lagr[sys])
    nnzRqr[sys]=np.count_nonzero(Rqr[sys])
    nnzGpercentual[sys]=nnzG[sys]*100/G[sys].size
    nnzLagpercentual[sys]=nnzLag[sys]*100/Lagr[sys].size
    nnzRqrpercentual[sys]=nnzRqr[sys]*100/Rqr[sys].size
#%%

i=0
for sys in lstsys:
    ax[i].spy(G[sys])
    ax[i].set_title("WLS-T "+sys)
    # ax[i].annotate("fnnz = {:.2f}%".format(nnzGpercentual[sys]),xycoords='subfigure fraction',xy=(0.0,0.0),xytext=(pos.x0-0.1,pos.y0-0.4), ha='center', va='center')
    ax[i].text(0.50,-0.07,"fnnz = {:.2f}%".format(nnzGpercentual[sys]),transform=ax[i].transAxes, ha='center', va='center')
    # ax[].annotate("fnnz = {:.2f}%".format(nnzLagpercentual),xycoords='subfigure fraction',xy=(0.43,0.28))
    # ax[i].grid(True,which="both")
    i=i+1
    ax[i].spy(Rqr[sys])
    
    ax[i].set_title("EnQR "+sys)
    ax[i].text(0.50,-0.07,"fnnz = {:.2f}%".format(nnzRqrpercentual[sys]),transform=ax[i].transAxes, ha='center', va='center')
    # ax[i].annotate("fnnz = {:.2f}%".format(nnzRqrpercentual[sys]),xycoords='subfigure fraction',xy=(0.0,0.0),xytext=(pos.x0,pos.y0))
    i=i+1

# plt.show()
plt.tight_layout()
fig.savefig("EsparsidadeSBSE"+sys+".pdf")

# %%
