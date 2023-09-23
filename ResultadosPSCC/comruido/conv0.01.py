#%%
#!/usr/bin/python
# -*- coding: UTF-8 -*-
import pandas as pd
import numpy as np
import numpy.linalg as liang
import scipy.sparse.linalg as sliang 




dsys={"IEEE14_rakp2009":"14 bus","IEEE118_rakp2009":"118 bus"}

dfResultados=pd.DataFrame()
for key,sys in dsys.items():
    dfres_x1_ComMed=pd.read_csv("resultados_conv_"+key+"-0.01ComMedidas.csv")
    dfres_x1_ComMed["ini"]=1
    dfres_x1_ComMed["Med"]=1
    dfres_x1_ComMed["sys"]=key
    dfres_x1_SeMMed=pd.read_csv("resultados_conv_"+key+"-0.01SemMedidas.csv")
    dfres_x1_SeMMed["ini"]=1
    dfres_x1_SeMMed["Med"]=0
    dfres_x1_SeMMed["sys"]=key
    dfResultados=pd.concat([dfResultados,dfres_x1_ComMed,dfres_x1_SeMMed])

#%%

dfResultados.to_csv("resultados_conv_completo_true.csv")
casos=list(set(dfResultados["caso"].tolist()))
#%%
dconvresults={}
dconvresults["caso"]=[]
dconvresults["Med"]=[]
dconvresults["ini"]=[]
dconvresults["sys"]=[]
dconvresults["convLM"]=[]
dconvresults["convGN"]=[]
dconvresults["convGNbc"]=[]

for key,sys in dsys.items():
    for caso in casos:
        for ini in [1]:
            for med in [0]:
                mask=(dfResultados["caso"]==caso)&(dfResultados["sys"]==key)& (dfResultados["ini"]==ini)& (dfResultados["Med"]==med) 
                nconvsLM=sum(dfResultados[mask]["convLM"].tolist())
                nconvsGN=sum(dfResultados[mask]["convGN"].tolist())
                nconvsGNbc=sum(dfResultados[mask]["convGNbc"].tolist())
                dconvresults["caso"].append(caso)
                dconvresults["sys"].append(key)
                dconvresults["ini"].append(ini)
                dconvresults["convLM"].append(nconvsLM)
                dconvresults["convGN"].append(nconvsGN)
                dconvresults["convGNbc"].append(nconvsGNbc)
                dconvresults["Med"].append(med)


dfResultadoFiltrados=pd.DataFrame(data=dconvresults)

#%%

dfResultadoFiltrados.to_csv("Resultado_conv_ruido_true.csv",index=None)
# %%
