#%%
#!/usr/bin/python
# -*- coding: UTF-8 -*-
import pandas as pd
import numpy as np
import numpy.linalg as liang
import scipy.sparse.linalg as sliang 
import matplotlib.pyplot as plt
import seaborn as sns


sns.color_palette()

dfDADOS=pd.read_csv("resultados_conv_completo.csv",index_col=0)
#%%

k=0.55
dsys={"IEEE14_rakp2009":"14 bus","IEEE118_rakp2009":"118 bus"}


fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(12*k,5*k))
#%%
i=0
AB=["a)","b)"]
for key, sys in dsys.items():
    
    ax[i].set_title(AB[i]+sys)
    boxplotLM=dfDADOS[(dfDADOS["convLM"]==1)&(dfDADOS["Med"]==0)&(dfDADOS["sys"]==key)]["nitsLM"].array
    boxplotGN=dfDADOS[(dfDADOS["convGN"]==1)&(dfDADOS["Med"]==0)&(dfDADOS["sys"]==key)]["nitsGN"].array
    boxplotGNbc=dfDADOS[(dfDADOS["convGNbc"]==1)&(dfDADOS["Med"]==0)&(dfDADOS["sys"]==key)]["nitsGNbc"].array
    dfLM=pd.DataFrame(data={"nº Its":boxplotLM})
    dfLM["method"]="LM"
    dfGN=pd.DataFrame(data={"nº Its":boxplotGN})
    dfGN["method"]="GN"
    dfGNbc=pd.DataFrame(data={"nº Its":boxplotGNbc})
    dfGNbc["method"]="GNbc"
    df=pd.concat([dfLM,dfGN,dfGNbc])
    ax[i].set_yticks(np.arange(0,11)*3)
    ax[i].grid(axis="y",zorder=3,ls=":")
    sns.boxplot(data=df,x="method",y="nº Its",zorder=2,notch=True,palette=["#099358","#D60D0D","#291C95"],showfliers=False,ax=ax[i])
    i=i+1
# %%
plt.tight_layout()
plt.savefig("Boxplot_semmeds.pdf")

# %%
