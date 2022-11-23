from classes import *
import numpy as np
import pandas as pd
from readfiles import *
import scipy.sparse.linalg as sliang 
import scipy.sparse as sparse 
from networkcalc import *
import numpy.linalg as liang

def calcCovRes(graph,dfDMED,ind_i):

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    W=create_W(z,flag_ones=0,prec_virtual=1e-5)
    calc_H_EE(z,var_t,var_v,graph,H)
    G=np.matmul(np.matmul(H.T,W),H)
    Ginv=liang.inv(G)
    S=np.matmul(np.matmul(H,Ginv),H.T)
    Winv=np.diag(1/np.diag(W))
    Cov=Winv-S
    return Cov


def renorm(graph,dfDMED,ind_i,cov):
    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    W=create_W(z,flag_ones=0,prec_virtual=1e-5)
    dz=np.zeros(len(z))
    calc_dz(z,graph,dz)
    Rn=np.abs(dz)/np.sqrt(np.diag(cov))
    bhat=np.sqrt(1/np.diag(W))*Rn*(1/np.sqrt(np.diag(cov)))
    zT=[]
    zde=[]
    zpara=[]
    for m in z:
        zT.append(m.type)
        zde.append(graph[m.k].bar.id)
        if m.type ==2 or m.type ==3:
            zpara.append(graph[m.m].bar.id)
        else:
            zpara.append(-1)
    d={"Tipo":zT,"de":zde,"para":zpara,"Res":dz,"Rn":Rn,"bhat":bhat}
    dfRes=pd.DataFrame(d)
    return dfRes
