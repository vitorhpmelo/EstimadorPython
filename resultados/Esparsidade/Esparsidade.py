
#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#%%
sys="IEEE118"

G=np.loadtxt(sys+"/G.csv",dtype=float,delimiter=",")
Lagr=np.loadtxt(sys+"/Lagra.csv",dtype=float,delimiter=",")
Rqr=np.loadtxt(sys+"/Rqr.csv",dtype=float,delimiter=",")


fig, (ax1,ax2,ax3)=plt.subplots(ncols=3)

nnzG=np.count_nonzero(G)
nnzLag=np.count_nonzero(Lagr)
nnzRqr=np.count_nonzero(Rqr)

nnzGpercentual=nnzG*100/G.size
nnzLagpercentual=nnzLag*100/Lagr.size
nnzRqrpercentual=nnzRqr*100/Rqr.size


ax1.spy(G)
ax1.set_title("G")
ax1.annotate("fnnz = {:.2f}%".format(nnzGpercentual),xycoords='subfigure fraction',xy=(0.11,0.28))
ax2.spy(Lagr)
ax2.set_title("Lagrangeano")
ax2.annotate("fnnz = {:.2f}%".format(nnzLagpercentual),xycoords='subfigure fraction',xy=(0.43,0.28))
ax3.spy(Rqr)
ax3.set_title("R-Fatoração QR")
ax3.annotate("fnnz = {:.2f}%".format(nnzRqrpercentual),xycoords='subfigure fraction',xy=(0.75,0.28))
plt.tight_layout()
fig.savefig("Esparsidade"+sys+".pdf")

# %%
