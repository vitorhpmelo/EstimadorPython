from classes import *
import numpy as np
import pandas as pd
from readfiles import *
import scipy.sparse.linalg as sliang 
import scipy.sparse as sparse 
from networkcalc import *
import numpy.linalg as liang
#file with the information of the libary

def NormalEQ(H,W,dz,prtG=0):
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    if(prtG==1):
        np.savetxt("G.csv",G,delimiter=",")
    A=sparse.csc_matrix(G)
    dx=sliang.spsolve(A,grad)
    return dx

def NormalEQ_CG(H,W,dz,prtG=0):
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    if(prtG==1):
        np.savetxt("G.csv",G,delimiter=",")
    A=sparse.csc_matrix(G)
    dx,flag=sliang.cg(A,grad)
    return dx



def NormalEQ_QR(H,W,dz,prtG=0):
    """
    
    """
    Wmei=np.diag(np.sqrt(np.diag(W)))
    H2=np.matmul(Wmei,H)
    [Q,R]=liang.qr(H2)
    b=np.matmul(np.matmul(Q.T,Wmei),dz)
    A=sparse.csr_matrix(R)
    dx=sliang.spsolve_triangular(A,b,lower=False)
    return dx




def SS_WLS(graph,dfDMED,ind_i,tol=1e-5,solver="QR"):
    """
    
    """
    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)

    Vinici(graph,flatStart=1)
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    dz=np.zeros(len(z))

    W=create_W(z,flag_ones=0,prec_virtual=1e-5)

    it=0
    while(it <10):
        calc_dz(z,graph,dz)
        calc_H_EE(z,var_t,var_v,graph,H)
        if(it==0 or it == 4):
            np.savetxt("H"+str(it)+".csv",H,delimiter=",")
        if solver=="Normal":
            dx=NormalEQ(H,W,dz)
        elif solver =="QR":
            dx=NormalEQ_QR(H,W,dz)
        elif solver == "cg":
            dx=NormalEQ_CG(H,W,dz)
        #dx=np.linalg.solve(G,grad)
        new_X(graph,var_t,var_v,dx)
        print("max dx {:e} ".format(np.amax(np.abs(dx))))
        if (np.amax(np.abs(dx))<tol):
            conv=1
            txt="Convergiu em {:d} iteracoes".format(it)
            print(txt)
            prt_state(graph)
            break
        it=it+1


def SS_WLS_lagrangian(graph,dfDMED,ind_i,tol=1e-5):
    Vinici(graph,flatStart=1)
    [z,c,var_t,var_v]=create_z_c_x_LGI(graph,dfDMED,ind_i)
    C=np.zeros((len(c),len(var_t)+len(var_v)))
    print("existem {:d} medidas virtuais".format(len(c)))
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    dz=np.zeros(len(z))
    cx=np.zeros(len(c))
    W=create_W(z,flag_ones=0,prec_virtual=1e-5)
    Zermat=np.zeros((C.shape[0],C.shape[0]))
    it=0
    while(it <10):
        calc_dz(z,graph,dz)
        calc_cx(c,graph,cx)
        calc_H_EE(z,var_t,var_v,graph,H)
        calc_H_EE(c,var_t,var_v,graph,C)
        grad=np.matmul(np.matmul(H.T,W),dz)
        G=np.matmul(np.matmul(H.T,W),H)
        b=np.concatenate((grad,-cx))
        M=np.concatenate((np.concatenate((G,C)),np.concatenate((C.T,Zermat))),axis=1)
        A=sparse.csc_matrix(M)
        dxl=sliang.spsolve(A,b)
        dx=dxl[:len(var_t)+len(var_v)]
        lamda=dxl[len(var_t)+len(var_v):]
        new_X(graph,var_t,var_v,dx)
        if (np.amax(np.abs(dx))<tol):
            conv=1
            txt="Convergiu em {:d} iteracoes".format(it)
            print(txt)
            prt_state(graph)
            break
        it=it+1

