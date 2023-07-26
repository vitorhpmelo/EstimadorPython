#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class bar():
    def __init__(self,id,type,contador):
        self.id=id
        self.type=type
        self.i=contador
        self.V=1
        self.teta=0
        self.Sbase=100
        self.Vbase=138
        self.Pg=0
        self.Qg=0
        self.Pd=0
        self.Qd=0
        self.Bs=0
        self.nloads=0
        self.nshunts=0
        self.ngds=0

class branch():
    def __init__(self,id,de,para,type,i):
        self.id=id
        self.de=de
        self.para=para
        self.type=type
        self.i=i
        self.x=-1
        self.r=-1
        self.ykm=0
        self.Y=np.zeros((2,2),dtype=complex)
        self.bsh=-1
        self.tap=-1
        self.limPA=-999
        self.flagLimP=0
    def cykm(self):
        self.ykm=1/complex(self.r,self.x)
    def twoPortCircuit(self):
        if self.type == 1:
            self.Y[0][0]=self.ykm+self.bsh
            self.Y[1][1]=self.ykm+self.bsh
            self.Y[1][0]=-self.ykm
            self.Y[0][1]=-self.ykm
        elif self.type ==2:
            self.Y[0][0]=((1/self.tap)**2)*self.ykm
            self.Y[1][1]=self.ykm
            self.Y[1][0]=-(1/self.tap)*self.ykm
            self.Y[0][1]=-(1/self.tap)*self.ykm
    def Pf(self,grafo,flagT):
        k=self.de
        m=self.para
        if flagT==0:
            P=(grafo[k].V**2)*np.real(self.Y[0][0]) + grafo[k].V* grafo[m].V*(\
            np.real(self.Y[0][1])*np.cos(grafo[k].teta-grafo[m].teta)\
            +np.imag(self.Y[0][1])*np.sin(grafo[k].teta-grafo[m].teta))
            return P
        elif flagT==1:
            P=(grafo[m].V**2)*np.real(self.Y[1][1]) + grafo[m].V* grafo[k].V*(\
            np.real(self.Y[1][0])*np.cos(grafo[m].teta-grafo[k].teta)\
            +np.imag(self.Y[0][1])*np.sin(grafo[m].teta-grafo[k].teta))
            return P
        else:
            return 0
    def Qf(self,grafo,flagT):
        k=self.de
        m=self.para
        if flagT==0:
            Qf=-(grafo[k].V**2)*np.imag(self.Y[0][0]) - grafo[k].V* grafo[m].V*(\
            np.imag(self.Y[0][1])*np.cos(grafo[k].teta-grafo[m].teta)\
            -np.real(self.Y[0][1])*np.sin(grafo[k].teta-grafo[m].teta))
            return Qf
        elif flagT==1:
            Qf=-(grafo[m].V**2)*np.imag(self.Y[1][1]) - grafo[m].V* grafo[k].V*(\
            np.imag(self.Y[1][0])*np.cos(grafo[m].teta-grafo[k].teta)\
            -np.real(self.Y[1][0])*np.sin(grafo[m].teta-grafo[k].teta))
            return Qf
        else:
            return 0                
    def dPfdt(self,grafo,flagT,var):
        k=self.de
        m=self.para
        Vk=grafo[k].V
        Vm=grafo[m].V
        tk=grafo[k].teta
        tm=grafo[m].teta
        if flagT==0:
            Bkm=np.imag(self.Y[0][1])
            Gkm=np.real(self.Y[0][1])
            if k==var: # dPkm/dtk
                return Vk*Vm*(Bkm*np.cos(tk-tm)-Gkm*np.sin(tk-tm))
            elif m==var:# dPkm/dtm
                return Vk*Vm*(-Bkm*np.cos(tk-tm)+Gkm*np.sin(tk-tm))
            else : 
                return 0
        elif flagT==1:
            Bmk=np.imag(self.Y[1][0])
            Gmk=np.real(self.Y[1][0])
            if k==var: # dPmk/dtk
                return Vk*Vm*(-Bmk*np.cos(tk-tm)-Gmk*np.sin(tk-tm))
            elif m==var: # dPmk/dtm
                return Vk*Vm*(Bmk*np.cos(tk-tm)+Gmk*np.sin(tk-tm))
            else : 
                return 0
    def dPfdV(self,grafo,flagT,var):
        k=self.de
        m=self.para
        Vk=grafo[k].V
        Vm=grafo[m].V
        tk=grafo[k].teta
        tm=grafo[m].teta  
        if flagT==0: #dPkm
            Gkk=np.real(self.Y[0][0])
            Bkm=np.imag(self.Y[0][1])
            Gkm=np.real(self.Y[0][1])
            if k==var: #dPkm/dVk             
                return 2*Gkk*Vk + Vm*(Bkm*np.sin(tk-tm)+Gkm*np.cos(tk-tm))
            elif m==var:
                return Vk*(Bkm*np.sin(tk-tm)+Gkm*np.cos(tk-tm))
            else:
                return 0    
        elif flagT==1:
            Gmm=np.real(self.Y[1][1])
            Bmk=np.imag(self.Y[1][0])
            Gmk=np.real(self.Y[1][0])
            if k==var: #dPkm/dVk  
                return Vm*(-Bmk*np.sin(tk-tm)+Gmk*np.cos(tk-tm))
            elif m==var:#dPkm/dVm
                return 2*Gmm*Vm + Vk*(-Bmk*np.sin(tk-tm)+Gmk*np.cos(tk-tm))
            else:
                return 0
    def dQfdt(self,grafo,flagT,var):
        k=self.de
        m=self.para
        Vk=grafo[k].V
        Vm=grafo[m].V
        tk=grafo[k].teta
        tm=grafo[m].teta  
        if flagT==0: #dQkm
            Bkm=np.imag(self.Y[0][1])
            Gkm=np.real(self.Y[0][1])
            if k==var: #dQkm/dtk
                return Vk*Vm*(Bkm*np.sin(tk-tm)+Gkm*np.cos(tk-tm))
            elif m==var: #dQkm/dtm
                return Vk*Vm*(-Bkm*np.sin(tk-tm)-Gkm*np.cos(tk-tm))
        elif flagT==1: #dQmk
            Bmk=np.imag(self.Y[1][0])
            Gmk=np.real(self.Y[1][0])
            if k==var: #dQmk/dtk
                return Vk*Vm*(Bmk*np.sin(tk-tm)-Gmk*np.cos(tk-tm))            
            elif m==var: #dQmk/dtm
                return Vk*Vm*(-Bmk*np.sin(tk-tm)+Gmk*np.cos(tk-tm))
            else:
                return 0
    def dQfdV(self,grafo,flagT,var):
        k=self.de
        m=self.para
        Vk=grafo[k].V
        Vm=grafo[m].V
        tk=grafo[k].teta
        tm=grafo[m].teta
        if flagT==0: #dQkm  
            Bkk=np.imag(self.Y[0][0])
            Bkm=np.imag(self.Y[0][1])
            Gkm=np.real(self.Y[0][1])
            if k==var: #dQkm/dvk
                return -2*Bkk*Vk+Vm*(-Bkm*np.cos(tk-tm)+Gkm*np.sin(tk-tm))
            elif m==var: #dQkm/dvm
                return Vk*(-Bkm*np.cos(tk-tm)+Gkm*np.sin(tk-tm))
            else:
                return 0
        elif flagT==1:
            Bmm=np.imag(self.Y[1][1])
            Bmk=np.imag(self.Y[1][0])
            Gmk=np.real(self.Y[1][0])
            if k==var: #dQmk/dvk
                return Vm*(-Bmk*np.cos(tk-tm)-Gmk*np.sin(tk-tm))
            elif m==var: #dQmk/dvm
                return -2*Bmm*Vm + Vk*(-Bmk*np.cos(tk-tm)-Gmk*np.sin(tk-tm))
            else:
                return 0
            

class branTCSC(branch):
    def __init__(self,id,de,para,type,i,a=1,xtcsc_ini=-1,Pfesp=0):
        super().__init__(id,de,para,type,i)
        self.a=a
        self.xtcsc_ini=xtcsc_ini
        self.xtcsc=xtcsc_ini
        self.Pfesp=Pfesp
    def AttY(self):
        self.Y[0][0]=complex(0,-1/self.xtcsc)
        self.Y[1][1]=complex(0,-1/self.xtcsc)
        self.Y[1][0]=complex(0,1/self.xtcsc)
        self.Y[0][1]=complex(0,1/self.xtcsc)
    def dPfdx(self,grafo,flagT):
        if flagT==0:
            k=self.de
            m=self.para
        elif flagT==1:
            k=self.para
            m=self.de
        return -grafo[k].V*grafo[m].V*((1/self.xtcsc)**2)*np.sin(grafo[k].teta-grafo[m].teta)
    def dQfdx(self,grafo,flagT):
        if flagT==0:
            k=self.de
            m=self.para
        elif flagT==1:
            k=self.para
            m=self.de
        return -(1/self.xtcsc**2)*((grafo[k].V**2)-grafo[k].V*grafo[m].V*np.cos(grafo[k].teta-grafo[m].teta))
    def dPfdB(self,grafo,flagT):
        if flagT==0:
            k=self.de
            m=self.para
        elif flagT==1:
            k=self.para
            m=self.de
        return grafo[k].V*grafo[m].V*np.sin(grafo[k].teta-grafo[m].teta)
    def dQfdB(self,grafo,flagT):
        if flagT==0:
            k=self.de
            m=self.para
        elif flagT==1:
            k=self.para
            m=self.de
        return (grafo[k].V**2)-grafo[k].V*grafo[m].V*np.cos(grafo[k].teta-grafo[m].teta)




class node_graph():
    def __init__(self,id,bar):
        self.V=1
        self.teta=0
        self.Bs=0
        self.adjk=dict()
        self.adjm=dict()
        self.ladjk=[]
        self.ladjm=[]
        self.SVC=None
        self.id=id
        self.bar=bar
        self.V=1
        self.teta=0
        self.FlagBS=0
        self.FlagTCSC=0
        self.FlagSVC=0
        self.FlagUPFC=0
        self.bFACTS_adjk=dict()
        self.bFACTS_adjm=dict()
        self.bUFPC_adjk=dict()
        self.bUFPC_adjm=dict()
    def P(self,graph):
        P=0
        if self.FlagSVC==1:
            P=P+self.SVC.Gk*self.V**2 
        for key,item in self.adjk.items():
            P=P+item.Pf(graph,0)
        for key,item in self.adjm.items():
            P=P+item.Pf(graph,1)
        for key,item in self.bUFPC_adjk.items():
            P=P+item.Pps(graph)
        for key,item in self.bUFPC_adjm.items():
            P=P+item.Psp(graph)
        if np.abs(P)<1e-12:
            P=0
        return P
    def Q(self,graph):
        if self.FlagBS==0:
            Q=0
        else:    
            Q=-self.Bs*self.V**2 
        if self.FlagSVC==1:
            Q=Q-self.SVC.Bk*self.V**2 
        for key,item in self.adjk.items():
            Q=Q+item.Qf(graph,0)
        for key,item in self.adjm.items():
            Q=Q+item.Qf(graph,1)
        for key,item in self.bUFPC_adjk.items():
            Q=Q+item.Qps(graph)
        for key,item in self.bUFPC_adjm.items():
            Q=Q+item.Qsp(graph)
        if np.abs(Q)<1e-12:
            Q=0
        return Q
    def dPdt(self,graph,bar):
        if self.i==bar:
            dPdt=0
            for key,item in self.adjk.items(): 
                dPdt=dPdt+item.dPfdt(graph,FlagT=0,var=bar)
            for key,item in self.adjm.items():
                dPdt=dPdt+item.dPfdt(graph,FlagT=0,var=bar)
            return dPdt
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdt(graph,0,bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdt(graph,1,bar)
        else:
            return 0
    def dPdV(self,graph,bar):
        dPdV=0
        if self.FlagSVC==1:
            dPdV=dPdV+2*self.SVC.Gk*self.V      
        if self.i==bar:
            for key,item in self.adjk.items(): 
                dPdV=dPdV+item.dPfdV(graph,FlagT=0,var=bar)
            for key,item in self.adjm.items():
                dPdV=dPdV+item.dPfdV(graph,FlagT=0,var=bar)
            return dPdV
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdV(graph,0,bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdV(graph,1,bar)
        else:
            return 0

    def dQdt(self,graph,bar):
        dQdt=0
        if self.i==bar:
            for key,item in self.adjk.items(): 
                dQdt=dQdt+item.dQdt(graph,FlagT=0,var=bar)
            for key,item in self.adjm.items():
                dQdt=dQdt+item.dQdt(graph,FlagT=0,var=bar)
            return dQdt
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dQdt(graph,0,bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return  self.adjk[str(self.i)+"-"+str(bar)].dQdt(graph,1,bar)
        else:
            return 0
    def dQdV(self,graph,bar): ## ENTRA AQUI A DERIVADA DO SVC
        if self.i==bar:
            if self.FlagBS==0:
                dQdV=0
            else:    
                dQdV=-2*self.Bs*self.V
            if self.FlagSVC==1:
                dQdV=dQdV-2*self.SVC.Bk*self.V      
            for key,item in self.adjk.items(): 
                dQdV=dQdV+item.dQdV(graph,FlagT=0,var=bar)
            for key,item in self.adjm.items():
                dQdV=dQdV+item.dQdV(graph,FlagT=0,var=bar)
            return dQdV
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return  self.adjk[str(self.i)+"-"+str(bar)].dQdV(graph,0,bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return  self.adjk[str(self.i)+"-"+str(bar)].dQdV(graph,1,bar)
        else:
            return  0  
        

class UPFC():
    def __init__(self,id,de,para,Vse_ini,t_se_ini,Vsh_ini,t_sh_ini,Psp,Qsp,Vp,Rse,Xse,Rsh,Xsh,Vse_max,Vse_min,Vsh_max,Vsh_min,mode):
        self.id=id #id do UPFC
        self.p=de # bus from, folliwing the graph order 
        self.s=para # bus to, folliwing the graph order 
        self.Vse_ini=Vse_ini # Vse_ini initialization for the series source
        self.Vse=Vse_ini # Vse series source voltage magnitude
        self.t_se_ini=t_se_ini*np.pi/180 # t_se_ini series source voltage phase angle initizalization value
        self.t_se=t_se_ini*np.pi/180 # t_se series source voltage phase angle
        self.Vsh_ini=Vsh_ini # Vsh_ini inirialization for the shunt source
        self.Vsh=Vsh_ini # Vsh shunt voltage source magnitude
        self.t_sh_ini=t_sh_ini*np.pi/180 # t_sh_ini shunt source voltage phase angle initizalization value
        self.t_sh=t_sh_ini*np.pi/180 # t_sh sgunt source voltage phase angle
        self.Psp_set=Psp # specified value for the active power flow over the UPFC 
        self.Qsp_set=Qsp # specified value for the reactve power flow over the UPFC
        self.Vp=Vp #specified value for the voltage magnitude in the terminal p
        self.Vs=1 # voltage magnitude at the "s" terminal
        self.Rse=Rse #resistence for the series source
        self.Xse=Xse #reactance for the series source
        self.Rsh=Rsh #resistence for the shunt source
        self.Xsh=Xsh #reactance for the shunt source
        self.Yse=1/complex(Rse,Xse)  #complex adimitance for the series source
        self.Ysh=1/complex(Rsh,Xsh) #complex adimitance for the shunt source
        self.gse=np.real(self.Yse)
        self.bse=np.imag(self.Yse)
        self.gsh=np.real(self.Ysh)
        self.bsh=np.imag(self.Ysh)
        self.Vse_max=Vse_max # voltage magintude superior limit fdor the series source
        self.Vse_min=Vse_min # voltage magintude inferior limit fdor the series source
        self.Vsh_max=Vsh_max # voltage magintude superior limit fdor the shunt source
        self.Vsh_min=Vsh_min # voltage magintude inferior limit fdor the shunt source
        self.mode=mode # mode 0 controles voltage at bus "p"/ 1 does note control voltage at bus "p"
    
    def Pps(self,graph):
        p=self.p
        s=self.s
        PartI=(graph[p].V**2)*(self.gse+self.gsh)
        PartII=-graph[p].V*graph[s].V*(self.gse*np.cos(graph[p].teta-graph[s].teta)+self.bse*np.sin(graph[p].teta-graph[s].teta))
        PartIII=-graph[p].V*self.Vse*(self.gse*np.cos(graph[p].teta-self.t_se)+self.bse*np.sin(graph[p].teta-self.t_se))
        PartIV=-graph[p].V*self.Vsh*(self.gsh*np.cos(graph[p].teta-self.t_sh)+self.bsh*np.sin(graph[p].teta-self.t_sh))
        return PartI+PartII+PartIII+PartIV
    def Qps(self,graph):
        p=self.p
        s=self.s
        PartI=(graph[p].V**2)*(self.bse+self.bsh)
        PartII=-graph[p].V*graph[s].V*(self.bse*np.cos(graph[p].teta-graph[s].teta)-self.gse*np.sin(graph[p].teta-graph[s].teta))
        PartIII=-graph[p].V*self.Vse*(self.bse*np.cos(graph[p].teta-self.t_se)-self.gse*np.sin(graph[p].teta-self.t_se))
        PartIV=-graph[p].V*self.Vsh*(self.bsh*np.cos(graph[p].teta-self.t_sh)-self.gsh*np.sin(graph[p].teta-self.t_sh))
        return -PartI-PartII-PartIII-PartIV
    def Psp(self,graph):
        p=self.p
        s=self.s
        PartIII=(graph[s].V**2)*(self.gse)
        PartI=-graph[s].V*graph[p].V*(self.gse*np.cos(graph[s].teta-graph[p].teta)+self.bse*np.sin(graph[s].teta-graph[p].teta))
        PartII=graph[s].V*self.Vse*(self.gse*np.cos(graph[s].teta-self.t_se)+self.bse*np.sin(graph[s].teta-self.t_se))
        return PartI+PartII+PartIII
    def Qsp(self,graph):
        p=self.p
        s=self.s
        PartI=-graph[s].V*graph[p].V*(self.bse*np.cos(graph[s].teta-graph[p].teta)-self.gse*np.sin(graph[s].teta-graph[p].teta))
        PartII=graph[s].V*self.Vse*(self.bse*np.cos(graph[s].teta-self.t_se)-self.gse*np.sin(graph[s].teta-self.t_se))
        PartIII=(graph[s].V**2)*(self.bse)
        return -PartI-PartII-PartIII
    def Pse(self,graph):
        """
        Active Power generated by series font
        """
        p=self.p
        s=self.s

        PartI=(self.Vse**2)*self.gse
        PartII=-self.Vse*graph[p].V*(self.gse*np.cos(self.t_se-graph[p].teta)+self.bse*np.sin(self.t_se-graph[p].teta))
        PartIII=self.Vse*graph[s].V*(self.gse*np.cos(self.t_se-graph[s].teta)+self.bse*np.sin(self.t_se-graph[s].teta))
        return PartI+PartII+PartIII
    
    def Qse(self,graph):
        """
        Reactive Power generated by series font
        """
        p=self.p
        s=self.s
        PartI=(self.Vse**2)*self.bse
        PartII=-self.Vse*graph[p].V*(self.bse*np.cos(self.t_se-graph[p].teta)-self.gse*np.sin(self.t_se-graph[p].teta))
        PartIII=self.Vse*graph[s].V*(self.bse*np.cos(self.t_se-graph[s].teta)-self.gse*np.sin(self.t_se-graph[s].teta))
        return -PartI-PartII-PartIII

    def Psh(self,graph):
        """
        Active Power generated by shunt font
        """
        p=self.p
        PartI=-(self.Vsh**2)*self.gsh
        PartII=self.Vsh*graph[p].V*(self.gsh*np.cos(self.t_sh-graph[p].teta)+self.bsh*np.sin(self.t_sh-graph[p].teta))
        return PartI+PartII
    
    def Qsh(self,graph):
        """
        Reactive Power generated by shunt font
        """
        p=self.p
        PartI=-(self.Vsh**2)*self.bsh
        PartII=self.Vsh*graph[p].V*(self.bsh*np.cos(self.t_sh-graph[p].teta)-self.gsh*np.sin(self.t_sh-graph[p].teta))
        return -PartI-PartII
    
    def dPpsdtp(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "p" (from) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vs*(bse*np.cos(tp - ts) - gse*np.sin(tp - ts))\
               -Vp*Vse*(bse*np.cos(tp - tse) - gse *np.sin(tp - tse)) \
            - Vp*Vsh*(bsh*np.cos(tp - tsh) - gsh*np.sin(tp - tsh))
    
    
    def dQpsdtp(self,graph):
        """
        Calcualte the derivative of the reactve power in respect to the "p" (from) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*Vs*((-bse)*np.sin(tp -ts) - gse*np.cos(tp - ts))\
              + Vp*Vse*((-bse)*np.sin(tp - tse) - gse*np.cos(tp - tse))\
              + Vp*Vsh*((-bsh)*np.sin(tp - tsh) - gsh*np.cos(tp - tsh))
    

    def dPpsdts(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "s" (to) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vs*(-bse*np.cos(tp - ts) + gse*np.sin(tp - ts))

    def dQpsdts(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (to) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*Vs*(bse*np.sin(tp - ts) + gse*np.cos(tp - ts))


    def dPpsdtse(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "se" (series voltage source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vse*(-bse*np.cos(tp - tse) + gse*np.sin(tp - tse))

    def dQpsdtse(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "se" (series voltage source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*Vse*(bse*np.sin(tp - tse) + gse*np.cos(tp - tse))

    def dPpsdtsh(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "sh" (shunt voltage source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vsh*(-bsh*np.cos(tp - tsh) + gsh*np.sin(tp - tsh))

    def dQpsdtsh(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "sh" (shunt voltage source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*Vsh*(bsh*np.sin(tp - tsh) + gsh*np.cos(tp - tsh))
    
    def dPpsdVp(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "p" (from terminal) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return 2*Vp*(gse + gsh) - Vs*(bse*np.sin(tp - ts) + gse*np.cos(tp - ts)) \
            - Vse*(bse*np.sin(tp- tse) + gse*np.cos(tp - tse)) \
            - Vsh*(bsh*np.sin(tp -tsh) + gsh*np.cos(tp -tsh))

    def dQpsdVp(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "p" (from terminal) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -2*Vp*(bse + bsh)+ \
        Vs*(bse*np.cos(tp - ts) - gse*np.sin(tp - ts)) + \
        Vse*(bse*np.cos(tp-tse) - gse*np.sin(tp- tse)) +\
        Vsh*(bsh*np.cos(tp - tsh) -gsh*np.sin(tp - tsh))


    def dPpsdVs(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (to terminal) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*(bse*np.sin(tp - ts) + gse*np.cos(tp - ts))

    def dQpsdVs(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (to bus) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*(bse*np.cos(tp - ts) - gse*np.sin(tp - ts))
    
    def dPpsdVse(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "se" (series source) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*(bse*np.sin(tp - tse) + gse*np.cos(tp - tse))

    def dQpsdVse(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (series source) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*(bse*np.cos(tp - tse) - gse*np.sin(tp - tse))

    def dPpsdVsh(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "sh" (shunt source) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*(bsh*np.sin(tp - tsh) + gsh*np.cos(tp - tsh))

    def dQpsdVsh(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "sh" (shunt voltage) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*(bsh*np.cos(tp - tsh) - gsh*np.sin(tp - tsh))

    #-------------------derivatives from s to p---------------------------------------------------

    def dPspdtp(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "p" (from) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vs*(-bse*np.cos(tp - ts) - gse*np.sin(tp - ts))
    
    
    def dQspdtp(self,graph):
        """
        Calcualte the derivative of the reactve power in respect to the "p" (from) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*Vs*(-bse*np.sin(tp - ts) + gse*np.cos(tp - ts))
    
    def dPspdts(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "s" (to) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vs*(bse*np.cos(tp- ts) + gse*np.sin(tp- ts)) \
            + Vs*Vse*(bse*np.cos(ts - tse) - gse*np.sin(ts - tse))

    def dQspdts(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (to) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*Vs*(bse*np.sin(tp - ts) - gse*np.cos(tp - ts)) -\
             Vs*Vse*(-bse*np.sin(ts - tse) - gse*np.cos(ts - tse))


    def dPspdtse(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "se" (series voltage source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vs*Vse*(-bse*np.cos(ts - tse) + gse*np.sin(ts - tse))

    def dQspdtse(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "se" (series voltage source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vs*Vse*(bse*np.sin(ts - tse) + gse*np.cos(ts - tse))

    def dPspdtsh(self,graph):
        """
        Calcualte the derivative of the active power in respect to the "sh" (shunt voltage source) voltage angle
        """    

        return 0

    def dQspdtsh(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "sh" (shunt voltage source) voltage angle
        """    

        return 0
    
    def dPspdVp(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "p" (from terminal) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vs*(-bse*np.sin(tp - ts) + gse*np.cos(tp - ts))

    def dQspdVp(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "p" (from terminal) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vs*(bse*np.cos(tp - ts) + gse*np.sin(tp - ts))


    def dPspdVs(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (to terminal) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*(-bse*np.sin(tp - ts) + gse*np.cos(tp - ts)) + 2*Vs*gse \
            + Vse*(bse*np.sin(ts - tse) + gse*np.cos(ts - tse))

    def dQspdVs(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (to bus) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vp*(bse*np.cos(tp - ts) + gse*np.sin(tp - ts)) - 2*Vs*bse \
            - Vse*(bse*np.cos(ts - tse) - gse*np.sin(ts - tse))
    
    def dPspdVse(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "se" (series source) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vs*(bse*np.sin(ts - tse) + gse*np.cos(ts - tse))

    def dQspdVse(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "s" (series source) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vs*(bse*np.cos(ts - tse) - gse*np.sin(ts - tse))

    def dPspdVsh(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "sh" (shunt source) voltage magnitude
        """    


        return 0

    def dQspdVsh(self,graph):
        """
        Calcualte the derivative of the reactive power in respect to the "sh" (shunt voltage) voltage magnitude
        """    

        return 0
    
    def dIgdtp(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "p" (from bus) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vse*(-bse*np.cos(tp - tse) - gse*np.sin(tp - tse)) -\
            Vp*Vsh*(-bsh*np.cos(tp - tsh)- gsh*np.sin(tp - tsh))

    
    def dIgdts(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "s" (to bus) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vs*Vse*(-bse*np.cos(ts - tse) - gse*np.sin(ts - tse))

    def dIgdtse(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "se" (series source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vse*(bse*np.cos(tp - tse) + gse*np.sin(tp - tse)) \
            + Vs*Vse*(bse*np.cos(ts - tse) + gse*np.sin(ts - tse))

    def dIgdtsh(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "sh" (shun source) voltage angle
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*Vsh*(bsh*np.cos(tp - tsh) + gsh*np.sin(tp - tsh))


    def dIgdVp(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "p" (from) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vse*(-bse*np.sin(tp - tse) + gse*np.cos(tp - tse)) \
            - Vsh*(-bsh*np.sin(tp - tsh) + gsh*np.cos(tp - tsh))

    def dIgdVs(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "s" (to bus) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return Vse*(-bse*np.sin(ts - tse) + gse*np.cos(ts - tse))

    def dIgdVse(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "se" (series source) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*(-bse*np.sin(tp - tse) + gse*np.cos(tp - tse)) \
            + Vs*(-bse*np.sin(ts - tse) + gse*np.cos(ts - tse)) + 2*Vse*gse

    def dIgdVsh(self,graph):
        """
        Calcualte the derivative of the power mismatch in respect to the "sh" (series source) voltage magnitude
        """    
        p=self.p
        s=self.s
        Vp=graph[p].V
        Vs=graph[s].V
        tp=graph[p].teta
        ts=graph[s].teta
        Vse=self.Vse
        Vsh=self.Vsh
        tse=self.t_se
        tsh=self.t_sh
        gse=self.gse
        bse=self.bse
        gsh=self.gsh
        bsh=self.bsh

        return -Vp*(-bsh*np.sin(tp - tsh) + gsh*np.cos(tp - tsh)) + 2*Vsh*gsh


class SVC():
    def __init__(self,id,bus,Rt,Xt,Bini,BMAX,BMIN,aini,amax,amin):
        self.id=id
        self.bus=bus
        self.Rt=Rt
        self.Xt=Xt
        self.Bini=Bini
        self.BSVC=Bini
        self.Xeq=(self.Xt-1/self.BSVC)
        self.Bk=-self.Xeq/(self.Xeq**2+self.Rt**2)
        self.Gk=self.Rt/(self.Xeq**2+self.Rt**2)
        self.BMAX=BMAX
        self.BMIN=BMIN
        self.aini=aini
        self.amax=amax
        self.amin=amin
    def attYk(self):
        self.Xeq=(self.Xt-1/self.BSVC)
        self.Bk=-self.Xeq/(self.Xeq**2+self.Rt**2)
        self.Gk=self.Rt/(self.Xeq**2+self.Rt**2)
    def dGkdBsvc(self):
        numerador=-2*self.Rt*self.BSVC*(self.Xt*self.BSVC-1)
        denominador=((self.Rt**2)*(self.BSVC**2)+(self.Xt*self.BSVC-1)**2)**2
        return numerador/denominador
    def dBkdBsvc(self):
        #parei aqui
        numerador=(self.Rt**2)*(self.BSVC**2)-(self.Xt**2)*(self.BSVC**2)+2*(self.Xt)*(self.BSVC)-1  
        denominador=((self.Rt**2)*(self.BSVC**2)+(self.Xt*self.BSVC-1)**2)**2
        return -numerador/denominador







class netinfo():
    def __init__(self,nbar,nram,nvar,nteta,nv) -> None:
        self.nbar=nbar
        self.nram=nram
        self.nvar=nvar
        self.nteta=nteta
        self.nv=nv

class meas():
    def __init__(self,k,m,type,val,prec) -> None:
        self.k=k
        self.m=m
        self.type=type
        self.val=val
        self.prec=prec
        self.sigma=np.abs(val)*prec/3
    def dz(self,graph):
        if self.type==0:
            return self.val-graph[self.k].P(graph) ## inserir fluxos do UPFC
        elif self.type==1:
            return self.val-graph[self.k].Q(graph) ## inserir fluxos do UPFC
        elif self.type==2:
            keyk=str(self.k)+"-"+str(self.m)
            keym=str(self.m)+"-"+str(self.k)
            if keyk in graph[self.k].adjk.keys():
                return self.val-graph[self.k].adjk[keyk].Pf(graph,0) # Se for UPFC branch
            elif keym in graph[self.k].adjm.keys():
                return self.val-graph[self.k].adjm[keym].Pf(graph,1) # Se for UPFC branch.
            elif keyk in graph[self.k].bUFPC_adjk.keys(): 
                return self.val-graph[self.k].bUFPC_adjk[keyk].Pps(graph)
            elif keym in graph[self.k].bUFPC_adjm.keys():
                return self.val-graph[self.k].bUFPC_adjm[keym].Psp(graph)
            else:
                print("medida de fluxo de potencia ativa com ramo não existente")
                exit(1)
        elif self.type==3:
            keyk=str(self.k)+"-"+str(self.m)
            keym=str(self.m)+"-"+str(self.k)
            if keyk in graph[self.k].adjk.keys():
                return self.val-graph[self.k].adjk[keyk].Qf(graph,0)
            elif keym in graph[self.k].adjm.keys():
                return self.val-graph[self.k].adjm[keym].Qf(graph,1)
            elif keyk in graph[self.k].bUFPC_adjk.keys(): 
                return self.val-graph[self.k].bUFPC_adjk[keyk].Qps(graph)
            elif keym in graph[self.k].bUFPC_adjm.keys():
                return self.val-graph[self.k].bUFPC_adjm[keym].Qsp(graph)
            else:
                print("medida de fluxo de potencia reativa com ramo não existente")
                exit(1)
        elif self.type==4:
            return self.val-graph[self.k].V
        else:
            print("Tipo de medida não existente")
            exit(1)
    def cx(self,graph):
        if self.type==0:
            return graph[self.k].P(graph)
        elif self.type==1:
            return graph[self.k].Q(graph)
        elif self.type==2:
            keyk=str(self.k)+"-"+str(self.m)
            keym=str(self.m)+"-"+str(self.k)
            if keyk in graph[self.k].adjk.keys():
                return graph[self.k].adjk[keyk].Pf(graph,0)
            elif keym in graph[self.k].adjm.keys():
                return graph[self.k].adjm[keym].Pf(graph,1)
            else:
                print("medida de fluxo de potencia ativa com ramo não existente")
                exit(1)
        elif self.type==3:
            keyk=str(self.k)+"-"+str(self.m)
            keym=str(self.m)+"-"+str(self.k)
            if keyk in graph[self.k].adjk.keys():
                return graph[self.k].adjk[keyk].Qf(graph,0)
            elif keym in graph[self.k].adjm.keys():
                return graph[self.k].adjm[keym].Qf(graph,1)
            else:
                print("medida de fluxo de potencia reativa com ramo não existente")
                exit(1)
        elif self.type==4:
            return graph[self.k].V
        else:
            print("Tipo de medida não existente")
            exit(1)




class state():
    def __init__(self,v,t) -> None:
        self.v=v
        self.t=t

