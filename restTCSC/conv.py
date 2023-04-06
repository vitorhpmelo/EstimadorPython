import pandas as pd
import matplotlib.pyplot as plt




dfcas1=pd.read_csv("convEECaso1cbkteDc.csv")
dfcas1B=pd.read_csv("convEEBCaso1cbkteDc.csv")

# dfcas2=pd.read_csv("convEEBCaso1cbktc3.csv")
# dfcas2B=pd.read_csv("convEEBCaso1sbktc2.csv")


#%%
# fig,ax = plt.subplots(nrows=2,ncols=2)


# ax[0][0].set_title("xini =-0.0100 dx")
# ax[0][0].set_xlabel("it")
# ax[0][0].set_ylabel("dx")
# ax[0][0].semilogy(range(len(dfcas1)),dfcas1["dx"],label="c\ backtracking ",color="r")
# ax[0][0].semilogy(range(len(dfcas1B)),dfcas1B["dx"],label="s\ backtracking",color="b")
# ax[0][0].grid()
# ax[0][0].legend()

# ax[0][1].set_title("xini =-0.0100 norm grad")
# ax[0][1].set_xlabel("it")
# ax[0][1].set_ylabel("norm grad")
# ax[0][1].semilogy(range(len(dfcas1)),dfcas1["dz"],label="c\ backtracking",color="r")
# ax[0][1].semilogy(range(len(dfcas1B)),dfcas1B["dz"],label="s\ backtracking",color="b")
# ax[0][1].grid()
# ax[0][1].legend()

# ax[1][0].set_title("xini =0.0100 dx")
# ax[1][0].set_xlabel("it")
# ax[1][0].set_ylabel("dx")
# ax[1][0].semilogy(range(len(dfcas2)),dfcas2["dx"],label="c\ backtracking",color="r")
# ax[1][0].semilogy(range(len(dfcas2B)),dfcas2B["dx"],label="s\ backtracking",color="b")
# ax[1][0].grid()
# ax[1][0].legend()

# ax[1][1].set_title("xini =0.0100 normgrad")
# ax[1][1].set_xlabel("it")
# ax[1][1].set_ylabel("norm grad")
# ax[1][1].semilogy(range(len(dfcas2)),dfcas2["dz"],label="c\ backtracking",color="r")
# ax[1][1].semilogy(range(len(dfcas2B)),dfcas2B["dz"],label="s\ backtracking",color="b")
# ax[1][1].grid()
# ax[1][1].legend()

# fig.tight_layout()
# plt.show()


# fig,ax = plt.subplots(nrows=2,ncols=2)


# ax[0][0].set_title("EE c/ B dx")
# ax[0][0].set_xlabel("it")
# ax[0][0].set_ylabel("dx")
# ax[0][0].semilogy(range(len(dfcas1)),dfcas1["dx"],label="-0.001",color="r")
# ax[0][0].semilogy(range(len(dfcas1B)),dfcas1B["dx"],label="0.001",color="b")
# ax[0][0].grid()
# ax[0][0].legend()

# ax[0][1].set_title("EE c/ B norm grad")
# ax[0][1].set_xlabel("it")
# ax[0][1].set_ylabel("norm grad")
# ax[0][1].semilogy(range(len(dfcas1)),dfcas1["dz"],label="-0.001",color="r")
# ax[0][1].semilogy(range(len(dfcas1B)),dfcas1B["dz"],label="0.001",color="b")
# ax[0][1].grid()
# ax[0][1].legend()

# ax[1][0].set_title("EE c/ B dx")
# ax[1][0].set_xlabel("it")
# ax[1][0].set_ylabel("dx")
# ax[1][0].semilogy(range(len(dfcas2)),dfcas2["dx"],label="-0.10",color="r")
# # ax[1][0].semilogy(range(len(dfcas2B)),dfcas2B["dx"],label="s\ backtracking",color="b")
# ax[1][0].grid()
# ax[1][0].legend()

# ax[1][1].set_title("EE c/ B normgrad")
# ax[1][1].set_xlabel("it")
# ax[1][1].set_ylabel("norm grad")
# ax[1][1].semilogy(range(len(dfcas2)),dfcas2["dz"],label="-0.10",color="r")
# # ax[1][1].semilogy(range(len(dfcas2B)),dfcas2B["dz"],label="s\ backtracking",color="b")
# ax[1][1].grid()
# ax[1][1].legend()

# fig.tight_layout()
# plt.show()

#%%

fig,ax = plt.subplots(nrows=1,ncols=2)


ax[0].set_title("dx")
ax[0].set_xlabel("it")
ax[0].set_ylabel("dx")
ax[0].semilogy(range(len(dfcas1)),dfcas1["dx"],label="EE c\ x ",color="r")
ax[0].semilogy(range(len(dfcas1B)),dfcas1B["dx"],label="EE c\ b",color="b")
ax[0].grid()
ax[0].legend()

ax[1].set_title("norm grad")
ax[1].set_xlabel("it")
ax[1].set_ylabel("norm grad")
ax[1].semilogy(range(len(dfcas1)),dfcas1["dz"],label="EE c\ x",color="r")
ax[1].semilogy(range(len(dfcas1B)),dfcas1B["dz"],label="EE c\ b",color="b")
ax[1].grid()
ax[1].legend()

plt.show()