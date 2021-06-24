#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 18:22:32 2021

@author: changkai
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from Regions import Regions
from Region_SIR import RGN_SIR
import itertools as it
import pandas as pd
import matplotlib.pyplot as plt


class SIR(Regions):
    
    # outflow coefficients
    # returns a dictionary with region names being the keys
    def OFC(self, dict_):
        mc = self.mobilityCoef(dict_)
        keyMat = self.keyMat()
        n = len(self.regions)
        ofc = {}
        for i in range(n):
            tempCS = 0
            tempCI = 0
            tempCR = 0
            for j in range(n):
                if i != j:
                    key = keyMat[i,j]
                    tempCS = tempCS + mc[key][3]
                    tempCI = tempCI + mc[key][4]
                    tempCR = tempCR + mc[key][5]

            name = self.regions[i].name
            ofc[name] = [tempCS, tempCI, tempCR]
        
        return ofc
    
     # build diagonal blocks
    def diagBlock(self, dict_, idx):
        ofc = self.OFC(dict_)
        name = self.regions[idx].name
        curOFC = ofc[name]
        curOFC = [-1 * item for item in curOFC] # take the negative value
        block = sp.diags(curOFC, 0, shape=[3,3], format = 'csr')
        
        return block
    
    # build a 3x3 block matrix
    def Block(self, dict_, row, col):
        mc = self.mobilityCoef(dict_)
        keyMat = self.keyMat()
        if row == col:
            block = self.diagBlock(dict_, row)
        else:
            key = keyMat[row, col]
            tempCoef = mc[key]
            block = sp.diags(tempCoef[0:3], 0, shape=[3,3], format='csr')
        return block
    
    # the linear matrix
    def AL(self, dict_):
        n = len(self.regions)
        row1 = self.Block(dict_, 0, 0)       
        
        for col in range(1,n):
            tempBlock = self.Block(dict_, 0, col)
            row1 = sp.hstack([row1, tempBlock])
        
        Al = row1
        for row in range(1,n):
            tempRow = self.Block(dict_, row, 0)
            for col in range(1,n):
                tempBlock = self.Block(dict_, row, col)
                tempRow = sp.hstack([tempRow, tempBlock])
            Al = sp.vstack([Al, tempRow])
        
        return Al
    
    # u is of shape (m,1); numpy array
    # return the dictionary containing the nonlinear terms of each region
    # in the order of S, I, R
    # def NLC(self, u):
    #     nlc = {}
    #     n = len(self.regions)
    #     ut = u
    #     m = ut.shape[0]
    #     ut = ut.reshape(1,m).squeeze()
    #     ut = ut.tolist()
    #     for i in range(n):
    #         name = self.regions[i].name
    #         nlc[name] = ut[3*i:(3*(i+1))]
    #     return nlc
    
    # return the nonlinear matrix (depending on u)
    def ANL(self, u):
        
        n = len(self.regions)
        B0 = sp.diags([0.], shape=[3, 3], format = 'csr') # zero block
        
        ut = u.reshape(1,3*n).squeeze()
        ut = ut.tolist()
    
        # the first row
        region1 = self.regions[0]
        S1, I1, R1 = ut[0:3]
        c1 = [-region1.beta*I1/region1.N, region1.beta*S1/region1.N-region1.gamma, 0]
        row1 = sp.diags(c1, 0, shape = [3,3], format = 'csr')
        row1[2,1] = region1.gamma
        for i in range(1,n):
            row1 = sp.hstack([row1, B0])

        # initialize the resultant matrix with the first row
        Anl = row1
        
        # loop over rows
        for i in range(1,n):
            # the first block of tempRow is a zero block
            tempRow = B0
            
            # stack each column
            for j in range(1,n):
                if i == j:
                    region = self.regions[i]
                    S, I, R = ut[3*i:3*(i+1)]
                    c = [-region.beta*I/region.N, region.beta*S/region.N-region.gamma, 0]
                    tempBlock = sp.diags(c, 0, shape = [3,3], format = 'csr')
                    tempBlock[2,1] = region.gamma
                else:
                    tempBlock = B0
                tempRow = sp.hstack([tempRow, tempBlock])
            
            # stack each row
            Anl = sp.vstack([Anl, tempRow])
            
        return Anl
    
    def Heun(self, dict_, dt, Tf, Ti=0):
        
        # define problem
        n = int((Tf - Ti)/dt)
        u0 = self.u0()
        Al = self.AL(dict_)
        Anl = self.ANL(u0)
        Au = Al + Anl
        u = u0
        # matrix that stores result
        U = u 
        for i in range(n):
            
            # predictor
            up = (dt*Au).dot(u)
        
            # corrector
            u = u + dt/2*(Au.dot(u) + Au.dot(up))
            
            # append result
            U = np.append(U, u, axis=1)
            
            # update ANL
            Anl = self.ANL(u)
            
            # update Au
            Au = Al + Anl
            
        return U
    
    def plot_Heun(self, dict_, index, dt, Tf, Ti=0, save = False):
        U = self.Heun(dict_, dt, Tf)
        n = int((Tf-Ti)/dt)+1
        tvals = np.linspace(Ti, Tf, n)
        fig, ax = plt.subplots()

        for j in range(3):
            if j == 0:
                Label = self.regions[index].name + '-S'
            elif j == 1:
                Label = self.regions[index].name + '-I'
            else:
                Label = self.regions[index].name + '-R'
            
            ax.plot(tvals, U[index*3+j,:], label = Label)
    
        ax.legend(loc = 'right')
        plt.show()
        if save:
            filename = self.regions[index].name + '.png'
            fig.savefig(filename)
    
    def plotAll_Heun(self, dict_, dt, Tf, Ti=0):
        U = self.Heun_SIR(dict_, dt, Tf)
        numRegion = len(self.regions)
        n = int((Tf-Ti)/dt)+1
        tvals = np.linspace(Ti, Tf, n)
        fig, ax = plt.subplots()
        for i in range(numRegion):
            for j in range(3):
                if j == 0:
                    Label = self.regions[i].name + '-S'
                elif j == 1:
                    Label = self.regions[i].name + '-I'
                else:
                    Label = self.regions[i].name + '-R'
            
            ax.plot(tvals, U[i,:], label = Label)
    
        ax.legend(loc = 'upper right')
        plt.show()
    
    
    
    
    
    
    
    