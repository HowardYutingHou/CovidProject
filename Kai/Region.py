#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 18:08:25 2021

@author: changkai
"""

class RGN(object):
    
    def __init__(self, Tf, Ti=0, name = 'unknown'):
        self.name = name
        self.Ti = Ti # initial time
        self.Tf = Tf # final time
        
    # build a key 
    def buildKey(self, region):
        return self.name + ' ' + region.name