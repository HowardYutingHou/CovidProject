#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 17:36:55 2021

@author: changkai
"""

import numpy as np 

class Regions(object):
    regions = []
    
    def __init__(self, *args):
        for region in args:
            self.regions.append(region)
            
    # return a matrix of keys
    def keyMat(self):
        n = len(self.regions)
        keyMat = np.zeros((n,n)).astype(str)
        for i in range(n):
            for j in range(n):
                keyMat[i,j] = self.regions[i].buildKey(self.regions[j])
        
        return keyMat
    
    # return a dictionary containing all available mobility coefficients
    # the order is mS_in, mI_in, mR_in, mS_out, mI_out, mR_out
    def mobilityCoef(self, dict_):
        keyMat = self.keyMat()
        n = len(self.regions)
        mc = {}
        for i in range(n):
            col = i+1
            for j in range(col, n):
                if keyMat[i,j] in dict_.keys():
                    curCoef = dict_[keyMat[i,j]]
                    mc[keyMat[i,j]] = curCoef
                    symCoef = curCoef[3:] + curCoef[0:3]
                    mc[keyMat[j,i]] = symCoef
                else:
                    if keyMat[j,i] not in dict_.keys():
                        raise RuntimeError('Dataset is not complete!\n'+
                                           'Need the mobility coefficients for '+
                                           keyMat[i,j] + ' pair!')
                    else:
                        symCoef = dict_[keyMat[j,i]]
                        mc[keyMat[j,i]] = symCoef
                        curCoef = symCoef[3:] + symCoef[0:3]
                        mc[keyMat[i,j]] = curCoef
        return mc
    
    # initial condition
    def u0(self):
        u0 = []
        for region in self.regions:
            u0 = u0 + region.y0
        n = len(u0)
        u0 = np.asarray(u0)
        u0 = u0.reshape((n,1))
        return u0
    
    
    
    
    
    
    
    
    
    
    
    
    