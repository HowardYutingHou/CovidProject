#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 18:13:52 2021
""
@author: changkai
"""
from Region import RGN

class RGN_SIR(RGN):
    
    def __init__(self, data, Tf, Ti=0):
        super().__init__(Tf, Ti=Ti, name = data['name'])
        self.y0 = [data['S'], data['I'], data['R']]
        self.N = sum(self.y0) # initial total population
        self.beta = data['beta']
        self.gamma = data['gamma']
