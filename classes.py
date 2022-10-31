#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class bar():
    def __init__(self,id,tipo,contador):
        self.id=id
        self.tipo=tipo
        self.i=contador
    V=1
    teta=0
    Sbase=100
    Vbase=138
    Pg=0
    Qg=0
    Pd=0
    Qd=0
    Bs=0
    nloads=0
    nshunts=0
    ngds=0

class branch():
    def __init__(self,id,de,para,tipo,i):
        self.id=id
        self.de=de
        self.para=para
        self.tipo=tipo
        self.i=i
    x=-1
    r=-1
    ykm=0
    Y=np.zeros((2,2),dtype=complex)
    bsh=-1
    tap=-1
    limPA=-999
    flagLimP=0
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
    adjk=dict()
    adjm=dict()
    V=1
    teta=0
    Bs=0
    def __init__(self,id,bar):
        self.id=id
        self.bar=bar    
        


class netinfo():
    def __init__(self,nbar,nram,nvar,nteta,nv) -> None:
        self.nbar=nbar
        self.nram=nram
        self.nvar=nvar
        self.nteta=nteta
        self.nv=nv