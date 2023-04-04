import pandas as pd
import matplotlib.pyplot as plt




dfcas1=pd.read_csv("convcaso1.csv")
dfcas1B=pd.read_csv("convcaso1_B.csv")

dfcas2=pd.read_csv("convcaso2.csv")
dfcas2B=pd.read_csv("convcaso2_B.csv")



fig,ax = plt.subplots(nrows=2,ncols=2)


ax[0][0].set_title("Caso 1 dx")
ax[0][0].set_xlabel("it")
ax[0][0].set_ylabel("dx")
ax[0][0].semilogy(range(len(dfcas1)),dfcas1["dx"],label="Using x",color="r")
ax[0][0].semilogy(range(len(dfcas1B)),dfcas1B["dx"],label="Using B",color="b")
ax[0][0].grid()
ax[0][0].legend()

ax[0][1].set_title("Caso 1 dz")
ax[0][1].set_xlabel("it")
ax[0][1].set_ylabel("dz")
ax[0][1].semilogy(range(len(dfcas1)),dfcas1["dz"],label="Using x",color="r")
ax[0][1].semilogy(range(len(dfcas1B)),dfcas1B["dz"],label="Using B",color="b")
ax[0][1].grid()
ax[0][1].legend()

ax[1][0].set_title("Caso 2 dx")
ax[1][0].set_xlabel("it")
ax[1][0].set_ylabel("dx")
ax[1][0].semilogy(range(len(dfcas2)),dfcas2["dx"],label="Using x",color="r")
ax[1][0].semilogy(range(len(dfcas2B)),dfcas2B["dx"],label="Using B",color="b")
ax[1][0].grid()
ax[1][0].legend()

ax[1][1].set_title("Caso 2 dz")
ax[1][1].set_xlabel("it")
ax[1][1].set_ylabel("dz")
ax[1][1].semilogy(range(len(dfcas2)),dfcas2["dz"],label="Using x",color="r")
ax[1][1].semilogy(range(len(dfcas2B)),dfcas2B["dz"],label="Using B",color="b")
ax[1][1].grid()
ax[1][1].legend()

fig.tight_layout()
plt.show()