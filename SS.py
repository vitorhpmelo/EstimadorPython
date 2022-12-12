from classes import *
import numpy as np
import pandas as pd
from readfiles import *
import scipy.sparse.linalg as sliang 
import scipy.sparse as sparse 
from networkcalc import *
import numpy.linalg as liang
import time as tm
#file with the information of the libary

def NormalEQ(H,W,dz,printcond=0,printmat=0):
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    if(printmat==1):
        np.savetxt("G.csv",G,delimiter=",",fmt="%.e")
    if(printcond==1):
        print("Ncond G(x) {:e}, Ncond H(x) {:e}".format(np.linalg.cond(G),np.linalg.cond(H)))
        with open("conds.csv","a") as f:
            f.write("{:e}\n".format(np.linalg.cond(G)))
    A=sparse.csc_matrix(G)
    try:
        dx=sliang.spsolve(A,grad)
    except:
        return -1
    return dx

def NormalEQ_CG(H,W,dz,printmat=0):
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    if(printmat==0):
        np.savetxt("G.csv",G,delimiter=",",fmt="%.e")
    A=sparse.csc_matrix(G)
    dx=np.zeros(G.shape[0])
    dx,flag=sliang.bicgstab(A,grad,dx,tol=1e-5,maxiter=20)
    print(dx)
    print("flag cg {:d}".format(flag))
    return dx



def NormalEQ_QR(H,W,dz,printcond=0,printmat=0):
    """
    
    """
    Wmei=np.diag(np.sqrt(np.diag(W)))
    H2=np.matmul(Wmei,H)
    if(printcond==1):
        print("Ncond WH(x) {:e}, Ncond H(x) {:e}".format(np.linalg.cond(H2),np.linalg.cond(H)))
        with open("conds.csv","a") as f:
            f.write("{:e}\n".format(np.linalg.cond(H2)))
    [Q,R]=liang.qr(H2)
    b=np.matmul(np.matmul(Q.T,Wmei),dz)
    A=sparse.csr_matrix(R)
    if printmat==1:
        np.savetxt('Rqr.csv',R,delimiter=",",fmt="%.e")
    dx=sliang.spsolve_triangular(A,b,lower=False)
    return dx




def SS_WLS(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-9,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,prinnormgrad=0):
    """
    
    """
    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)

    Vinici(graph,flatStart=1)
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    dz=np.zeros(len(z))

    W=create_W(z,flag_ones=0,prec_virtual=prec_virtual)
    backtracking=1
    it=0
    tit=[]
    ts=tm.time()

    f=open('convIEEE'+str(len(graph))+solver+str(int(-np.log10(prec_virtual)))+".csv","w")

    while(it <9):
        t1=tm.time()
        calc_dz(z,graph,dz)
        calc_H_EE(z,var_t,var_v,graph,H)
        grad=np.matmul(np.matmul(H.T,W),dz)
        print("Norma do gradiente {:e}".format(liang.norm(grad)))
        if it==0 and prinnormgrad==1:
            norminicial=liang.norm(grad)
        if(it==0 or it == 4):
            np.savetxt("H"+str(it)+".csv",H,delimiter=",")
        if solver=="Normal":
            if it==0:
                dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
            else:
                dx=NormalEQ(H,W,dz,printcond=printcond,printmat=0)
                if len(dx)==1:
                    break ##matrix singular
        elif solver =="QR":
            if it==0:
                dx=NormalEQ_QR(H,W,dz,printcond=printcond,printmat=printmat)
            else:
                dx=NormalEQ_QR(H,W,dz,printcond=printcond,printmat=0)
        elif solver == "cg":
            dx=NormalEQ_CG(H,W,dz,printmat=printmat)
        #dx=np.linalg.solve(G,grad)
        
        new_X(graph,var_t,var_v,dx)
        #fbacktracking(graph,dx,z,var_t,var_v,H,dz,W)

        t2=tm.time()
        tit.append(t2-t1)
        print("max dx {:e} ".format(np.amax(np.abs(dx))))
        if prinnormgrad==1:
            calc_dz(z,graph,dz)
            calc_H_EE(z,var_t,var_v,graph,H)
            grad=np.matmul(np.matmul(H.T,W),dz)
            f.write("{:d},{:.3e},{:.3}\n".format(it,liang.norm(grad),np.amax(np.abs(dx))))
            if liang.norm(grad)/norminicial < tol2:
                txt="Convergiu em {:d} iteracoes".format(it)
                print(txt)
                prt_state(graph)
                break
        if prinnormgrad!=1:
            if (np.amax(np.abs(dx))<tol):
                conv=1
                txt="Convergiu em {:d} iteracoes".format(it)
                print(txt)
                # prt_state(graph)
                break
        it=it+1
    tf=tm.time()
    f.close()
    return (tf-ts),tit



def SS_WLS_lagrangian(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-9,printcond=0,printmat=0,printnormgrad=0):
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
    tit=[]
    ts=tm.time()
    f=open('convIEEE'+str(len(graph))+"lagran.csv","w")

    while(it <10):
        t1=tm.time()
        calc_dz(z,graph,dz)
        calc_cx(c,graph,cx)
        calc_H_EE(z,var_t,var_v,graph,H)
        calc_H_EE(c,var_t,var_v,graph,C)
        grad=np.matmul(np.matmul(H.T,W),dz)
        if printnormgrad==1:
            print("Norma do gradiente {:e}".format(liang.norm(np.abs(grad))))
            if it==0:
                norminicial=liang.norm(grad)
        G=np.matmul(np.matmul(H.T,W),H)
        b=np.concatenate((grad,-cx))
        M=np.concatenate((np.concatenate((G,C)),np.concatenate((C.T,Zermat))),axis=1)
        if printcond==1:
            print("Ncond HLa(x) {:e}, Ncond H(x) {:e}".format(np.linalg.cond(M),np.linalg.cond(H)))
            with open("conds.csv","a") as file:
                file.write("{:e}\n".format(np.linalg.cond(M)))
        if printmat==1:
            np.savetxt("Lagra.csv",M,delimiter=",",fmt="%.e")
        A=sparse.csc_matrix(M)
        dxl=sliang.spsolve(A,b)
        dx=dxl[:len(var_t)+len(var_v)]
        lamda=dxl[len(var_t)+len(var_v):]
        new_X(graph,var_t,var_v,dx)
        t2=tm.time()
        tit.append(t2-t1)
        if printnormgrad==1:
            calc_dz(z,graph,dz)
            calc_H_EE(z,var_t,var_v,graph,H)
            grad=np.matmul(np.matmul(H.T,W),dz)
            f.write("{:d},{:.3e},{:.3e}\n".format(it,liang.norm(grad),np.amax(np.abs(dx))))
            if liang.norm(grad)/norminicial < tol2:
                txt="Convergiu em {:d} iteracoes".format(it)
                print(liang.norm(grad)/norminicial)
                print(txt)
                prt_state(graph)
                break
        else:    
            if (np.amax(np.abs(dx))<tol):
                conv=1
                txt="Convergiu em {:d} iteracoes".format(it)
                print(txt)
                prt_state(graph)
                break
        it=it+1
    tf=tm.time()
    f.close()
    return tf-ts,tit


def fbacktracking(graph,dx,z,var_t,var_v,H,dz,W):
    calc_dz(z,graph,dz)
    calc_H_EE(z,var_t,var_v,graph,H)
    grad=np.matmul(np.matmul(H.T,W),dz)
    objtk=np.matmul(dz,np.matmul(W,dz))
    new_X(graph,var_t,var_v,dx)
    objtnk=np.matmul(dz,np.matmul(W,dz))
    it=0
    while(objtnk < objtk+1e-4*np.dot(grad,dx)):
        new_X(graph,var_t,var_v,-dx)
        dx=dx/2
        calc_dz(z,graph,dz)
        new_X(graph,var_t,var_v,dx)
        objtnk=np.matmul(dz,np.matmul(W,dz))
        if it>10:
            print("backtrackin falhou")
        break
    
def get_state(graph):
    v=[]
    teta=[]
    for no in graph:
        v.append(no.V)
        teta.append(no.teta)
    v=np.array(v)
    teta=np.array(teta)
    state(v,teta)
    return state(v,teta)


def SS_WLS_clean(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-9,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,prinnormgrad=0):
    """
    Withouth printing options
    """
    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)

    Vinici(graph,flatStart=1)
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    dz=np.zeros(len(z))

    W=create_W(z,flag_ones=0,prec_virtual=prec_virtual)
    backtracking=1
    it=0
    tit=[]
    ts=tm.time()

    while(it <9):
        t1=tm.time()
        calc_dz(z,graph,dz)
        calc_H_EE(z,var_t,var_v,graph,H)
        grad=np.matmul(np.matmul(H.T,W),dz)
        if it==0 and prinnormgrad==1:
            norminicial=liang.norm(grad)
        if solver=="Normal":
            if it==0:
                dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
            else:
                dx=NormalEQ(H,W,dz,printcond=printcond,printmat=0)
                if len(dx)==1:
                    break ##matrix singular
        elif solver =="QR":
            if it==0:
                dx=NormalEQ_QR(H,W,dz,printcond=printcond,printmat=printmat)
            else:
                dx=NormalEQ_QR(H,W,dz,printcond=printcond,printmat=0)
        elif solver == "cg":
            dx=NormalEQ_CG(H,W,dz,printmat=printmat)
        #dx=np.linalg.solve(G,grad)
        
        new_X(graph,var_t,var_v,dx)
        #fbacktracking(graph,dx,z,var_t,var_v,H,dz,W)

        t2=tm.time()
        tit.append(t2-t1)
        if prinnormgrad==1:
            calc_dz(z,graph,dz)
            calc_H_EE(z,var_t,var_v,graph,H)
            grad=np.matmul(np.matmul(H.T,W),dz)
            if liang.norm(grad)/norminicial < tol2:
                conv=1
                break
        if prinnormgrad!=1:
            if (np.amax(np.abs(dx))<tol):
                conv=1
                break
        it=it+1
    tf=tm.time()
    print("convergência {:d}".format(conv))
    return (tf-ts),tit,conv,it



def SS_WLS_lagrangian_clean(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-9,printcond=0,printmat=0,printnormgrad=0):
    """
    Withouth printing options
    """
    Vinici(graph,flatStart=1)
    [z,c,var_t,var_v]=create_z_c_x_LGI(graph,dfDMED,ind_i)
    C=np.zeros((len(c),len(var_t)+len(var_v)))
    H=np.zeros((len(z),len(var_t)+len(var_v)))
    dz=np.zeros(len(z))
    cx=np.zeros(len(c))
    W=create_W(z,flag_ones=0,prec_virtual=1e-5)
    Zermat=np.zeros((C.shape[0],C.shape[0]))
    it=0
    tit=[]
    ts=tm.time()

    while(it <10):
        t1=tm.time()
        calc_dz(z,graph,dz)
        calc_cx(c,graph,cx)
        calc_H_EE(z,var_t,var_v,graph,H)
        calc_H_EE(c,var_t,var_v,graph,C)
        grad=np.matmul(np.matmul(H.T,W),dz)
        if printnormgrad==1:
            if it==0:
                norminicial=liang.norm(grad)
        G=np.matmul(np.matmul(H.T,W),H)
        b=np.concatenate((grad,-cx))
        M=np.concatenate((np.concatenate((G,C)),np.concatenate((C.T,Zermat))),axis=1)
        A=sparse.csc_matrix(M)
        dxl=sliang.spsolve(A,b)
        dx=dxl[:len(var_t)+len(var_v)]
        lamda=dxl[len(var_t)+len(var_v):]
        new_X(graph,var_t,var_v,dx)
        t2=tm.time()
        tit.append(t2-t1)
        if printnormgrad==1:
            calc_dz(z,graph,dz)
            calc_H_EE(z,var_t,var_v,graph,H)
            grad=np.matmul(np.matmul(H.T,W),dz)
            if liang.norm(grad)/norminicial < tol2:
                # txt="Convergiu em {:d} iteracoes".format(it)
                # print(liang.norm(grad)/norminicial)
                # print(txt)
                # prt_state(graph)
                conv=1
                break
        else:    
            if (np.amax(np.abs(dx))<tol):
                conv=1
                # txt="Convergiu em {:d} iteracoes".format(it)
                # print(txt)
                # prt_state(graph)
                break
        it=it+1
    tf=tm.time()
    print("convergência {:d}".format(conv))
    return tf-ts,tit,conv,it
