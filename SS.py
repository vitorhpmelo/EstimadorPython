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
        np.savetxt("G.csv",G,delimiter=",",fmt="%.15e")
    if(printcond==1):
        # print("Ncond G(x) {:e}, Ncond H(x) {:e}".format(np.linalg.cond(G),np.linalg.cond(H)))
        with open("conds.csv","a") as f:
            f.write("{:e},{:e}\n".format(np.linalg.cond(G),np.linalg.cond(H)))
    A=sparse.csc_matrix(G)
    try:
        dx=sliang.spsolve(A,grad)
    except:
        return -1
    return dx


def NormalEQ_lev(H,W,damp,dz,printcond=0,printmat=0,D=1):
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    # I=np.diag(np.diag(G))
    I=np.eye(len(G))

    I=np.eye(len(G))*D
    A=np.add(G,damp*I)

    if(printmat==1):
        np.savetxt("G.csv",G,delimiter=",",fmt="%.15e")
    if(printcond==1):
        # print("Ncond G(x) {:e}, Ncond H(x) {:e}".format(np.linalg.cond(G),np.linalg.cond(H)))
        with open("conds.csv","a") as f:
            f.write("{:e},{:e}\n".format(np.linalg.cond(G),np.linalg.cond(H)))
    A=sparse.csc_matrix(A)
    try:
        dx=sliang.spsolve(A,grad)
    except:
        return -1
    return dx

def NormalEQ_lev_scale(H,W,damp,dz,printcond=0,printmat=0):
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    I=np.eye(len(G))*(np.max(np.diag(G)))
    
    A=np.add(G,damp*I)

    if(printmat==1):
        np.savetxt("G.csv",G,delimiter=",",fmt="%.15e")
    if(printcond==1):
        # print("Ncond G(x) {:e}, Ncond H(x) {:e}".format(np.linalg.cond(G),np.linalg.cond(H)))
        with open("conds.csv","a") as f:
            f.write("{:e},{:e}\n".format(np.linalg.cond(G),np.linalg.cond(H)))
    A=sparse.csc_matrix(A)
    try:
        dx=sliang.spsolve(A,grad)
    except:
        return -1
    return dx


def NormalEQ_CG(H,W,dz,printmat=0):
    grad=np.matmul(np.matmul(H.T,W),dz)
    G=np.matmul(np.matmul(H.T,W),H)
    if(printmat==0):
        np.savetxt("G.csv",G,delimiter=",",fmt="%.15e")
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
        np.savetxt('Rqr.csv',R,delimiter=",",fmt="%.15e")
    dx=sliang.spsolve_triangular(A,b,lower=False)
    return dx




def SS_WLS(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-9,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,prinnormgrad=0):
    """
    Função que executa o estimador de estado WLS com diferentes sovers.
    @solver == "QR utiliza a fatoração QR
    @solver == "Normal utiliza a equação normal
    @solver == "cg" utiliza gradientes conjugados (fase de testes)
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

    while(it <10):
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
                prt_state(graph)
                break
        it=it+1
    tf=tm.time()
    f.close()
    return (tf-ts),tit



def SS_WLS_lagrangian(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-9,printcond=0,printmat=0,printnormgrad=0):
    """
    Função que perfoma a estimação com igualdades lagrangianas
    """
    
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
            np.savetxt("Lagra.csv",M,delimiter=",",fmt="%.15e")
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
    
def get_state(graph,sample="ref"):
    d={}
    d["tipo"]=[]
    d["de"]=[]
    d["val"]=[]
    d["sample"]=[]
    for no in graph:
        #v
        d["tipo"].append("v")
        d["de"].append(no.id)
        d["val"].append(no.V)
        d["sample"].append(sample)
        #teta
        d["tipo"].append("teta")
        d["de"].append(no.id)
        d["val"].append(no.teta)
        d["sample"].append(sample)
    dfAns=pd.DataFrame(data=d)
    return dfAns


def get_state_TCSC(ramTCSC):
    x={}
    for key,ram in ramTCSC.items():
        x[key]=float(ram.xtcsc)
    return x

def get_state_FACTS(TCSC={},svc={},UPFC={},sample="ref"):
    d={}
    d["tipo"]=[]
    d["de"]=[]
    d["val"]=[]
    d["sample"]=[]
    for key,ram in TCSC.items():
        d["tipo"].append("x_tcsc")
        d["de"].append(key)
        d["val"].append(ram.xtcsc)
        d["sample"].append(sample)
    for key,s in svc.items():
        d["tipo"].append("B_svc")
        d["de"].append(key)
        d["val"].append(s.BSVC)
        d["sample"].append(sample)
    for key,u in UPFC.items():
        ##Vsh
        d["tipo"].append("UPFC_Vsh")
        d["de"].append(key)
        d["val"].append(u.Vsh)
        d["sample"].append(sample)
        ##tsh
        d["tipo"].append("UPFC_tsh")
        d["de"].append(key)
        d["val"].append(u.t_sh)
        d["sample"].append(sample)
        #vse
        d["tipo"].append("UPFC_Vse")
        d["de"].append(key)
        d["val"].append(u.Vse)
        d["sample"].append(sample)
        #tse
        d["tipo"].append("UPFC_tse")
        d["de"].append(key)
        d["val"].append(u.t_se)
        d["sample"].append(sample)

    dfANS=pd.DataFrame(data=d)
    return dfANS

def SS_WLS_clean(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-9,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,prinnormgrad=0):
    """
    Withouth printing options for computing time
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
    conv=0
    while(it <20):
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
    Withouth printing options for computing time
    """

    conv=0
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

    while(it <20):
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




def SS_WLS_FACTS(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V-0.1

    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=2
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <30):
        a=1
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)
        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)
        grad=np.matmul(np.matmul(H.T,W),b)
        # dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
        dx=NormalEQ_QR(H,W,b,printcond=printcond,printmat=printmat)
        Jxk=np.matmul(np.matmul(b,W),b)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
            calc_dz(z,graph,dz)
            calc_cUPFC(graph,var_UPFC,c_upfc)
            b=np.append(dz,c_upfc)
            Jxn=np.matmul(np.matmul(b,W),b)
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-a*dx)
                new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-a*dx)
                a=a/2
                it2=it2+1
        print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            print(liang.norm(grad)/norminicial)
            print(txt)
            prt_state(graph)
            prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break


        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv_A.csv', index=False)

    return conv

def SS_WLS_FACTS_noBC(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printgrad=1,printres=1,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V-0.01



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=2
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <30):
        a=1
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)
        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)
        grad=np.matmul(np.matmul(H.T,W),b)
<<<<<<< HEAD
        try:
            dx=NormalEQ_QR(H,W,b,printcond=printcond,printmat=printmat)
=======
        try: 
            dx=NormalEQ(H,W,b,printcond=printcond,printmat=printmat)
>>>>>>> refs/remotes/origin/TCSC+SVC+UPFC
        except:
            conv=0
            it=30
            break
<<<<<<< HEAD
            
=======

>>>>>>> refs/remotes/origin/TCSC+SVC+UPFC
        Jxk=np.matmul(np.matmul(b,W),b)
        if it==0:
            norminicial=liang.norm(grad)
        new_X(graph,var_t,var_v,a*dx)
        new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
        new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
        new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        b=np.append(dz,c_upfc)
        Jxn=np.matmul(np.matmul(b,W),b)

        if printgrad==True:
            print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)
        lstdx.append(maxdx)
<<<<<<< HEAD

=======
        lstdz.append(gradredux)
>>>>>>> refs/remotes/origin/TCSC+SVC+UPFC
        if maxdx>1e3:
            conv=0
            it=30
            break
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            if printres==True:
                print(liang.norm(grad)/norminicial)
                print(txt)
                prt_state(graph)
                prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break

        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        dfits = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        dfits.to_csv('conv_GN.csv', index=False)
    elif pirntits==2:
        iterdict={"dx":lstdx,"dz":lstdz}
        dfits = pd.DataFrame(iterdict)
    else:
        dfits=[]
    return conv,it,dfits

def SS_WLS_FACTS_withBC(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printres=1,printgrad=1,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[0])
            graph[m].V=graph[m].V-0.01


    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=5
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <30):
        a=1
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)
        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)
        grad=np.matmul(np.matmul(H.T,W),b)
        try: 
            dx=NormalEQ(H,W,b,printcond=printcond,printmat=printmat)
        except:
            conv=0
            it=30
            break

        # dx=NormalEQ_QR(H,W,b,printcond=printcond,printmat=printmat)
        Jxk=np.matmul(np.matmul(b,W),b)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
            # new_X_EE_UPFC_lim(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,dx,a)
            calc_dz(z,graph,dz)
            calc_cUPFC(graph,var_UPFC,c_upfc)
            b=np.append(dz,c_upfc)
            Jxn=np.matmul(np.matmul(b,W),b)
            it2=it2+1
            if it2==itmax:
                break
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-a*dx)
                new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-a*dx)
                # new_X_EE_UPFC_lim(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-dx,a)
                a=a/2
        if printgrad==True:   
            print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)

        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if maxdx>1e3:
            conv=0
            it=30
            break
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            if printres==True:
                print(liang.norm(grad)/norminicial)
                print(txt)
                prt_state(graph)
                prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break


        it=it+1


    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        dfits = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        dfits.to_csv('conv_GNbc.csv', index=False)
    elif pirntits==2:
        iterdict={"dx":lstdx,"dz":lstdz}
        dfits = pd.DataFrame(iterdict)
    else:
        dfits=[]

    return conv,it,dfits

def SS_WLS_FACTS_grad(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[0])
            graph[m].V=graph[m].V+0.1



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=5
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
    
    while(it <30):
        a=1
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)
        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)
        grad=np.matmul(np.matmul(H.T,W),b)
        if it<5:
            dx=grad/np.linalg.norm(grad)
        else:
            dx=NormalEQ_QR(H,W,b,printcond=printcond,printmat=printmat)
        Jxk=np.matmul(np.matmul(b,W),b)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
            calc_dz(z,graph,dz)
            calc_cUPFC(graph,var_UPFC,c_upfc)
            b=np.append(dz,c_upfc)
            Jxn=np.matmul(np.matmul(b,W),b)
            print("funçao objetivo k {} | k+1 : {}".format(Jxk,Jxn))
            it2=it2+1
            if it2==itmax:
                break
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-a*dx)
                new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-a*dx)
                if it<5:
                    a=a/20
                else:
                    a=a/2
                
        prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
        print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            print(liang.norm(grad)/norminicial)
            print(txt)
            prt_state(graph)
            prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break


        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv_A.csv', index=False)

    return it


def SS_WLS_FACTS_LM(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printgrad=1,printres=1,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices LevenberMerquard

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V-0.1



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=2
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <30):
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)

        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)

        grad=np.matmul(np.matmul(H.T,W),b)
        Jxk=np.matmul(np.matmul(b,W),b)

        if it==0:
            norminicial=liang.norm(grad)
            Jxkin=Jxk
            G=np.matmul(np.matmul(H.T,W),H)
            D=liang.norm(np.diag(G))*0.0001

            
        damp=calc_damp_leven_mod_2(grad/norminicial,it+1)
        dx=NormalEQ_lev(H,W,damp,b,D=D)

        new_X(graph,var_t,var_v,dx)
        new_X_TCSC(graph,len(var_t)+len(var_v),var_x,dx)
        new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,dx)
        new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,dx)
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        b=np.append(dz,c_upfc)
        if printgrad==True:
            print("grad {:e}, dx {:e}".format( liang.norm(grad)/norminicial,liang.norm(dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            if printres==True:
                print(liang.norm(grad)/norminicial)
                print(txt)
                prt_state(graph)
                prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break

        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv_A.csv', index=False)

    return conv,it



def SS_WLS_FACTS_LM_BC(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printgrad=1,printres=1,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices LevenberMerquard

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V-0.01



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=1
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <35):
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)

        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)

        grad=np.matmul(np.matmul(H.T,W),b)
        Jxk=np.matmul(np.matmul(b,W),b)

        if it==0:
            norminicial=liang.norm(grad)
            G=np.matmul(np.matmul(H.T,W),H)
            D=liang.norm(np.diag(G))*1e-6

            
        damp=calc_damp_leven_mod_2(grad/norminicial,it+1)
        dx=NormalEQ_lev(H,W,damp,b,D=D)
        a=1
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
            calc_dz(z,graph,dz)
            calc_cUPFC(graph,var_UPFC,c_upfc)
            b=np.append(dz,c_upfc)
            Jxn=np.matmul(np.matmul(b,W),b)
            it2=it2+1
            if it2==itmax:
                break
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-a*dx)
                new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-a*dx)
                a=a/2
                
        if printgrad==True:
            print("grad {:e}, dx {:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            if printres==True:
                print(liang.norm(grad)/norminicial)
                print(txt)
                prt_state(graph)
                prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break

        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        dfits = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        dfits.to_csv('conv_LM.csv', index=False)
    elif pirntits==2:
        iterdict={"dx":lstdx,"dz":lstdz}
        dfits = pd.DataFrame(iterdict)
    else:
        dfits=[]

    return conv,it,dfits

def SS_WLS_FACTS_LM_3(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices LevenberMerquard

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V-0.1



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=2
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    LD=3
    LI=2
    count=0
    beta=2
    gama=3
    
    v=beta

    
    while(it <100):
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)

        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)

        grad_xk=np.matmul(np.matmul(H.T,W),b)
        Jxk=np.matmul(np.matmul(b,W),b)

        if it==0:
            norminicial=liang.norm(grad_xk)
            Jxkin=Jxk
            G=np.matmul(np.matmul(H.T,W),H)
            damp=np.max(np.diag(G))*0.001
        # gradnormal=grad/norminicial
        # Jxknorma=Jxk/Jxkin
        # damp=calc_damp_leven_mod(gradnormal)
        

        dx=NormalEQ_lev(H,W,damp,b)
        
        m0=cal_model_quad(grad_xk,Jxk,np.zeros(len(dx)),H,damp,np.eye(len(dx)),W)
        mk=cal_model_quad(grad_xk,Jxk,dx,H,damp,np.eye(len(dx)),W)
        model=cal_model_2(grad_xk,dx,damp)
        new_X(graph,var_t,var_v,dx)
        new_X_TCSC(graph,len(var_t)+len(var_v),var_x,dx)
        new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,dx)
        new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,dx)
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        b=np.append(dz,c_upfc)
        Jxk_novo=np.matmul(np.matmul(b,W),b)
        
        pk=(Jxk-Jxk_novo)/(model)
        # if pk<1/4:
        #     count=0
        #     damp=damp*LI
        #     if pk<0:
        #         new_X(graph,var_t,var_v,-dx)
        #         new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-dx)
        #         new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-dx)
        #         new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-dx)
        # elif pk>3/4:
        #     damp=damp/LD
        #     count=count+1
        #     if count>2:
        #         damp=damp/LD
        if pk>0:
            A=1/gama
            B=1-(beta-1)*(2*pk-1)**3
            damp=damp*np.max([A,B])
            v=beta
        else:
            damp=damp*v
            v=2*v
            new_X(graph,var_t,var_v,-dx)
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-dx)






        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)

        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)



        grad=np.matmul(np.matmul(H.T,W),b)
        print("grad {:e}, dx {:e}".format( liang.norm(grad)/norminicial,liang.norm(dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            print(liang.norm(grad)/norminicial)
            print(txt)
            prt_state(graph)
            prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break

        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv_A.csv', index=False)

    return it






def calc_damp_leven(grad,Jxk):
    return (np.linalg.norm(grad)**2)/Jxk

def calc_damp_leven_mod(grad):
    return (np.linalg.norm(grad))

def calc_damp_leven_mod_2(grad,it):
    return 2*np.linalg.norm(grad)/(10*(it**4))

def cal_model_quad(grad,Jx,dx,H,damp,D,W):

    return Jx+np.dot(dx,grad)+0.5*np.dot(dx,np.matmul(np.add(np.matmul(np.matmul(H.T,W),H),damp*D),dx))
    
def cal_model_2(grad,dx,damp):

    return -0.5*np.dot(dx,np.add(damp*dx,-grad))
    


def SS_WLS_FACTS_withBC_limalphavarfacts(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[0])
            graph[m].V=graph[m].V+0.1



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=5
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <30):
        a=1
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)
        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)
        grad=np.matmul(np.matmul(H.T,W),b)
        # dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
        dx=NormalEQ_QR(H,W,b,printcond=printcond,printmat=printmat)
        Jxk=np.matmul(np.matmul(b,W),b)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            if 0.1<dx_TCSC_max(graph,len(var_t)+len(var_v),var_x,dx):
                X_TCSC_its(graph,len(var_t)+len(var_v),var_x,dx)    
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
            calc_dz(z,graph,dz)
            calc_cUPFC(graph,var_UPFC,c_upfc)
            b=np.append(dz,c_upfc)
            Jxn=np.matmul(np.matmul(b,W),b)
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-a*dx)
                new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-a*dx)
                a=a/2
                it2=it2+1
        print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            print(liang.norm(grad)/norminicial)
            print(txt)
            prt_state(graph)
            prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break


        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv_A.csv', index=False)

    return it


def SS_WLS_FACTS_withBC_itvarfacts(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[0])
            graph[m].V=graph[m].V+0.1



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=3
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <30):
        a=1
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)
        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)
        grad=np.matmul(np.matmul(H.T,W),b)
        # dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
        dx=NormalEQ_QR(H,W,b,printcond=printcond,printmat=printmat)
        Jxk=np.matmul(np.matmul(b,W),b)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            if it<2:
                X_TCSC_its(graph,len(var_t)+len(var_v),var_x,dx)
            if it==1:
                reini_X_TCSC(graph,var_x,z)
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
            calc_dz(z,graph,dz)
            calc_cUPFC(graph,var_UPFC,c_upfc)
            b=np.append(dz,c_upfc)
            Jxn=np.matmul(np.matmul(b,W),b)
            it2=it2+1
            if it2==itmax:
                break
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-a*dx)
                new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-a*dx)
                a=a/2
        print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            print(liang.norm(grad)/norminicial)
            print(txt)
            prt_state(graph)
            prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break


        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv_A.csv', index=False)

    return it




def SS_WLS_FACTS_withBC_limvarfacts(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    var_svc=create_x_SVC(graph)
    [var_UPFC,c_upfc]=create_c_x_UPFC(graph)
    #create var UPFC

    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[0])
            graph[m].V=graph[m].V+0.1



    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    HSVC=np.zeros((len(z),len(var_svc)))
    UPFC=np.zeros((len(z),4*len(var_UPFC)))
    n_teta=len(var_t)
    n_v=len(var_v)
    n_TCSC=len(var_x)
    n_SVC=len(var_svc)
    n_UPFC=len(var_UPFC)
    nvar=n_teta+n_v+n_TCSC+n_SVC+4*n_UPFC
    dz=np.zeros(len(z))
    W=create_W(z+list(c_upfc),flag_ones=0,prec_virtual=prec_virtual) #expandir W para caber as c_FACTS
    
    C_UPFC=np.zeros((len(c_upfc),nvar))

    it=0
    it2=0
    itmax=5
    lstdx=[]
    lstdz=[]
    lstc_upfc=[]
    
    while(it <30):
        a=1
        calc_dz(z,graph,dz)
        calc_cUPFC(graph,var_UPFC,c_upfc)
        calc_H_EE(z,var_t,var_v,graph,Htrad) 
        calc_H_EE_TCSC(z,var_x,graph,HTCSC) 
        calc_H_EE_SVC(z,var_svc,graph,HSVC) 
        calc_H_EE_UPFC(z,var_UPFC,graph,UPFC)
        calc_C_EE_UPFC(var_t,var_v,var_x,var_svc,var_UPFC,graph,C_UPFC)
        
        Hx=np.concatenate((Htrad,HTCSC,HSVC,UPFC),axis=1)
        H=np.concatenate((Hx,C_UPFC),axis=0)
        b=np.append(dz,c_upfc)
        grad=np.matmul(np.matmul(H.T,W),b)
        # dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
        dx=NormalEQ_QR(H,W,b,printcond=printcond,printmat=printmat)
        Jxk=np.matmul(np.matmul(b,W),b)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            if 0.1<dx_TCSC_max(graph,len(var_t)+len(var_v),var_x,dx):
                X_TCSC_its(graph,len(var_t)+len(var_v),var_x,dx)    
            new_X_TCSC(graph,len(var_t)+len(var_v),var_x,a*dx)
            new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,a*dx)
            new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,a*dx)
            
            calc_dz(z,graph,dz)
            calc_cUPFC(graph,var_UPFC,c_upfc)
            b=np.append(dz,c_upfc)
            Jxn=np.matmul(np.matmul(b,W),b)
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                checklim_X_TCSC(graph,var_x)
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                new_X_SVC(graph,len(var_t)+len(var_v)+len(var_x),var_svc,-a*dx)
                new_X_EE_UPFC(graph,len(var_t)+len(var_v)+len(var_x)+len(var_svc),var_UPFC,-a*dx)
                a=a/2
                it2=it2+1
        print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if gradredux <tol2 and maxdx<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            upfc_angle(graph)
            print(liang.norm(grad)/norminicial)
            print(txt)
            prt_state(graph)
            prt_state_FACTS(graph,var_x,var_svc,var_UPFC)
            conv=1
            break


        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)

        # Save the DataFrame to a CSV file
        df.to_csv('conv_A.csv', index=False)

    return it

def SS_WLS_FACTS_clean(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,printmat=0,pirntits=0,prinnormgrad=0,flatstart=-1):
    
    '''
    WLS state estimator with FACTS devices (only TCSC implemented yet)

    @param graph with the informations of the network
    @param prt param indicating if it is printing everyting or not
    @param tol tolerance for the dx atualization of the variables
    @param tol2 tolerance for the gradiente reduction
    @param solver only gain matrix implemented yet
    @param prec_virtual standard deviation of virtual measurements
    @param printcond flag for calculating and printing condition number
    @param printmat flag for calculating and printing the matrix for calculationg the descend direction
    @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
    the flat start, 1 it uses the DBAR
    '''

    tin=tm.time()
    tits=[]
    conv=0
    c1=1e-4 #constant for backintracking
    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V+0.1
    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    dz=np.zeros(len(z))
    W=create_W(z,flag_ones=0,prec_virtual=prec_virtual)
    it=0
    it2=0
    itmax=10
    condlst=[]
    while(it <20):
        t0=tm.time()
        a=1
        calc_dz(z,graph,dz)
        calc_H_EE(z,var_t,var_v,graph,Htrad)
        calc_H_EE_TCSC(z,var_x,graph,HTCSC)
        H=np.concatenate((Htrad,HTCSC),axis=1)
        grad=np.matmul(np.matmul(H.T,W),dz)
        dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
        # dx=NormalEQ_QR(H,W,dz)
        Jxk=np.matmul(np.matmul(dz,W),dz)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            new_X_TCSCC(graph,len(var_t)+len(var_v),var_x,a*dx)
            calc_dz(z,graph,dz)
            Jxn=np.matmul(np.matmul(dz,W),dz)
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSCC(graph,len(var_t)+len(var_v),var_x,-a*dx)
                a=a/2
                it2=it2+1
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)

        if gradredux <tol2 and maxdx<tol:
            conv=1
            t1=tm.time()
            tits.append(t1-t0)
            break

        it=it+1
        t1=tm.time()
        tits.append(t1-t0)

    tf=tm.time()
    return conv,it,tits,tf-tin



def SS_WLS_FACTS_2(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,pirntits=0,printmat=0,prinnormgrad=0,flatstart=-1):
   
    '''

        WLS state estimator with FACTS devices (only TCSC implemented yet), but now using B as the derivative

        @param graph with the informations of the network
        @param prt param indicating if it is printing everyting or not
        @param tol tolerance for the dx atualization of the variables
        @param tol2 tolerance for the gradiente reduction
        @param solver only gain matrix implemented yet
        @param prec_virtual standard deviation of virtual measurements
        @param printcond flag for calculating and printing condition number
        @param printmat flag for calculating and printing the matrix for calculationg the descend direction
        @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
        the flat start, 1 it uses the DBAR

    '''
    
    c1=1e-4

    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V+0.1
    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    dz=np.zeros(len(z))
    W=create_W(z,flag_ones=0,prec_virtual=prec_virtual)
    it=0
    it2=0
    itmax=10

    lstdx=[]
    lstdz=[]
    while(it <20):
        a=1
        calc_dz(z,graph,dz)
        calc_H_EE(z,var_t,var_v,graph,Htrad)
        calc_H_EE_TCSC_B(z,var_x,graph,HTCSC)
        H=np.concatenate((Htrad,HTCSC),axis=1)
        grad=np.matmul(np.matmul(H.T,W),dz)
        dx=NormalEQ(H,W,dz,printcond=printcond,printmat=printmat)
        Jxk=np.matmul(np.matmul(dz,W),dz)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            new_X_TCSCC_B(graph,len(var_t)+len(var_v),var_x,a*dx)
            calc_dz(z,graph,dz)
            Jxn=np.matmul(np.matmul(dz,W),dz)
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSCC_B(graph,len(var_t)+len(var_v),var_x,-a*dx)
                a=a/2
            it2=it2+1
        print("{:e},{:e}".format( liang.norm(grad)/norminicial,liang.norm(a*dx)))    
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)
        lstdx.append(maxdx)
        lstdz.append(gradredux)
        if liang.norm(grad)/norminicial<tol2 and liang.norm(a*dx)<tol:
            txt="Convergiu em {:d} iteracoes".format(it)
            print(liang.norm(grad)/norminicial)
            print(txt)
            prt_state(graph)
            break


        it=it+1

    if pirntits==1:
        iterdict={"dx":lstdx,"dz":lstdz}
        df = pd.DataFrame(iterdict)
        # Save the DataFrame to a CSV file
        df.to_csv('conv_B.csv', index=False)

def SS_WLS_FACTS_2_clean(graph,dfDMED,ind_i,tol=1e-7,tol2=1e-7,solver="QR",prec_virtual=1e-5,printcond=0,pirntits=0,printmat=0,prinnormgrad=0,flatstart=-1):
   
    '''

        WLS state estimator with FACTS devices (only TCSC implemented yet), but now using B as the derivative

        @param graph with the informations of the network
        @param prt param indicating if it is printing everyting or not
        @param tol tolerance for the dx atualization of the variables
        @param tol2 tolerance for the gradiente reduction
        @param solver only gain matrix implemented yet
        @param prec_virtual standard deviation of virtual measurements
        @param printcond flag for calculating and printing condition number
        @param printmat flag for calculating and printing the matrix for calculationg the descend direction
        @param flat start, initialization of the state variables, if -1 uses the DC state estimator to intialize the angles and the X, 0 it ujses
        the flat start, 1 it uses the DBAR

    '''
    
    tin=tm.time()
    tits=[]
    conv=0
    c1=1e-4

    FACTSini(graph)

    Vinici(graph,flatStart=flatstart,dfDMED=dfDMED,ind_i=ind_i)

    [z,var_t,var_v]=create_z_x(graph,dfDMED,ind_i)
    var_x=create_x_TCSC(graph)
    if flatstart==2:
        for key in var_x.keys():
            key=key.split("-")
            m=int(key[1])
            graph[m].V=graph[m].V+0.1
    Htrad=np.zeros((len(z),len(var_t)+len(var_v)))
    HTCSC=np.zeros((len(z),len(var_x)))
    dz=np.zeros(len(z))
    W=create_W(z,flag_ones=0,prec_virtual=prec_virtual)
    it=0
    it2=0
    itmax=5


    while(it <30):
        t0=tm.time()
        a=1
        calc_dz(z,graph,dz)
        calc_H_EE(z,var_t,var_v,graph,Htrad)
        calc_H_EE_TCSC_B(z,var_x,graph,HTCSC)
        H=np.concatenate((Htrad,HTCSC),axis=1)
        grad=np.matmul(np.matmul(H.T,W),dz)
        dx=NormalEQ(H,W,dz,printcond=0,printmat=0)
        Jxk=np.matmul(np.matmul(dz,W),dz)
        if it==0:
            norminicial=liang.norm(grad)
        it2=0
        while it2<itmax:
            new_X(graph,var_t,var_v,a*dx)
            new_X_TCSCC_B(graph,len(var_t)+len(var_v),var_x,a*dx)
            calc_dz(z,graph,dz)
            Jxn=np.matmul(np.matmul(dz,W),dz)
            if Jxn < Jxk + c1*a*np.dot(grad,dx):
                break
            else:
                new_X(graph,var_t,var_v,-a*dx)
                new_X_TCSCC_B(graph,len(var_t)+len(var_v),var_x,-a*dx)
                a=a/2
            it2=it2+1
        gradredux=liang.norm(grad)/norminicial
        maxdx= liang.norm(a*dx)

        if gradredux<tol2 and maxdx<tol:
            conv=1
            t1=tm.time()
            tits.append(t1-t0)
            break
        it=it+1
        t1=tm.time()
        tits.append(t1-t0)

    tf=tm.time()
    return conv,it,tits,tf-tin



