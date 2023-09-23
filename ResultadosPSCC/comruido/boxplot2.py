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
dmeas={0:"FACTS not measured",1:"FACTS Measured"}

dsys={"IEEE118_rakp2009":"118 bus"}
sys="IEEE118_rakp2009"


fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(12*k,5*k))
#%%
i=0
AB=["a) ","b) "]

for key, meas in dmeas.items(): 
    ax[i].set_title(AB[i]+meas)
    boxplotLM=dfDADOS[(dfDADOS["convLM"]==1)&(dfDADOS["Med"]==key)&(dfDADOS["sys"]==sys)]["nitsLM"].array
    boxplotGN=dfDADOS[(dfDADOS["convGN"]==1)&(dfDADOS["Med"]==key)&(dfDADOS["sys"]==sys)]["nitsGN"].array
    boxplotGNbc=dfDADOS[(dfDADOS["convGNbc"]==1)&(dfDADOS["Med"]==key)&(dfDADOS["sys"]==sys)]["nitsGNbc"].array
    dfLM=pd.DataFrame(data={"nº Its":boxplotLM})
    dfLM["method"]="LM"
    dfGN=pd.DataFrame(data={"nº Its":boxplotGN})
    dfGN["method"]="GN"
    dfGNbc=pd.DataFrame(data={"nº Its":boxplotGNbc})
    dfGNbc["method"]="GNbc"
    df=pd.concat([dfLM,dfGN,dfGNbc])
    ax[i].set_yticks(np.arange(0,11)*3)
    ax[i].grid(axis="y",zorder=3,ls=":")
    sns.boxplot(data=df,x="method",y="nº Its",zorder=2,notch=True,palette=["#099358","#ff6961","#add8e6"],showfliers=False,ax=ax[i])
    i=i+1
# %%
plt.tight_layout()
plt.savefig("Boxplot_opt2.pdf")

# %%
