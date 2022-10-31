#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class bar():
    def __init__(self,id,tipo,contador):
        self.id=id
        self.tipo=tipo
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
    def __init__(self,id,de,para,tipo,i):
        self.id=id
        self.de=de
        self.para=para
        self.tipo=tipo
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
        if self.tipo == 1:
            self.Y[0][0]=self.ykm+self.bsh
            self.Y[1][1]=self.ykm+self.bsh
            self.Y[1][0]=-self.ykm
            self.Y[0][1]=-self.ykm
        elif self.tipo ==2:
            self.Y[0][0]=((1/self.tap)**2)*self.ykm
            self.Y[1][1]=self.ykm
            self.Y[1][0]=-(1/self.tap)*self.ykm
            self.Y[0][1]=-(1/self.tap)*self.ykm

class node_graph():

    def __init__(self,id,bar):
        self.V=1
        self.teta=0
        self.Bs=0
        self.adjk=dict()
        self.adjm=dict()
        self.id=id
        self.bar=bar    
        


class netinfo():
    def __init__(self,nbar,nram,nvar,nteta,nv) -> None:
        self.nbar=nbar
        self.nram=nram
        self.nvar=nvar
        self.nteta=nteta
        self.nv=nv