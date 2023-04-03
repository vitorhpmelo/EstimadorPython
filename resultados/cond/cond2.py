#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


color1="#391DC2"
color2="#BF1E17"
color3="#188C6C"
color4="#FBC519"

#%%

k=0.40
plt.title("k por Iteração Sistema de 118 barras")
plt.subplots(figsize=(k*16,k*5))
Norm=np.loadtxt("conitNorm.csv",dtype=float)
Lagra=np.loadtxt("conitLagra.csv",dtype=float)
plt.semilogy(range(len(Norm)),Norm,label="WLS-T",marker="d",color=color1)
plt.semilogy(range(len(Lagra)),Lagra,label="IgL",marker="o",color=color2)
plt.xlabel("Iteração")
plt.ylabel("k")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("condit.pdf")

# %%
