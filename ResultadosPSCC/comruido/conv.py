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
    dfres_x1_ComMed=pd.read_csv("resultados_conv_"+key+"x1ComMedidas.csv")
    dfres_x1_ComMed["ini"]=1
    dfres_x1_ComMed["Med"]=1
    dfres_x1_ComMed["sys"]=key
    dfres_x2_ComMed=pd.read_csv("resultados_conv_"+key+"x2ComMedidas.csv")
    dfres_x2_ComMed["ini"]=2
    dfres_x2_ComMed["Med"]=1
    dfres_x2_ComMed["sys"]=key
    dfres_x1_SeMMed=pd.read_csv("resultados_conv_"+key+"x1SemMedidas.csv")
    dfres_x1_SeMMed["ini"]=1
    dfres_x1_SeMMed["Med"]=0
    dfres_x1_SeMMed["sys"]=key
    dfres_x2_SeMMed=pd.read_csv("resultados_conv_"+key+"x2SemMedidas.csv")
    dfres_x2_SeMMed["ini"]=2
    dfres_x2_SeMMed["Med"]=0
    dfres_x2_SeMMed["sys"]=key
    dfResultados=pd.concat([dfResultados,dfres_x1_ComMed,dfres_x1_SeMMed,dfres_x2_ComMed,dfres_x2_SeMMed])

#%%

dfResultados.to_csv("resultados_conv_completo.csv")
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
        for ini in [1,2]:
            for med in [1,0]:
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

dfResultadoFiltrados.to_csv("Resultado_conv_ruido.csv",index=None)