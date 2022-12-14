#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%


plt.title("Número de Condicionamento por Iteração Sistema de 118 barras")
Norm=np.loadtxt("conitNorm.csv",dtype=float)
Lagra=np.loadtxt("conitLagra.csv",dtype=float)
plt.semilogy(range(len(Norm)),Norm,label="Condicionamento G",marker="d",color="r")
plt.semilogy(range(len(Lagra)),Lagra,label="Condicionamento Lagrangiano",marker="o",color="b")
plt.xlabel("Iteração")
plt.ylabel("Nº de Condicionamento")
plt.grid()
plt.legend()
plt.savefig("condit.pdf")

# %%
