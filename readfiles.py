from classes import *
import pandas as pd
import numpy as np

def read_files(sys):
    try:
        dfDBAR=pd.read_csv(sys+"/DBAR.csv",header=None,dtype={0:np.int64,1:np.int64})
        dfDBAR.columns=["id","type","V","teta","Pg","Qg","Pd","Qd","Bs"]
    except:
        print("Error while reading DBAR file")
        exit(1)

    try:
        dfDBRAN=pd.read_csv(sys+"/DBRAN.csv",header=None,dtype={0:np.int64,1:np.int64,2:np.int64,3:np.int64})
        dfDBRAN.columns=["id","type","de","para","r","x","bsh","tap"]
    except:
        print("Error while reading DBRAN file")
        exit(1)
    try:
        dfDMED=pd.read_csv(sys+"/DMED.csv",header=None)
    except:
        print("There is no DMED")
        dfDMED=[]
    return dfDBAR,dfDBRAN,dfDMED
