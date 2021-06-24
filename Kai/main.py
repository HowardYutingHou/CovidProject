#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 16:32:49 2021

@author: changkai
"""
#%% Import and Adjust Data
from SIR import SIR
from Region_SIR import RGN_SIR
import pandas as pd

# mobility data
df_m = pd.read_csv('/Users/changkai/Desktop/mathResearch/AV/COVID/TestingData/fake_inter.csv')
                   # names = ['Region1', 'Region2', 'S_in', 'I_in', 'R_in', 'S_out', 'I_out', 'R_out'])
# region data
df_r = pd.read_csv('/Users/changkai/Desktop/mathResearch/AV/COVID/TestingData/fake_data.csv')
                   # names = ['Region', 'S', 'I', 'R', 'beta', 'gamma', 'Name'])

# dict_ = df_m.set_index('Region1').T.to_dict('list')
df_m['Between'] = df_m['Region1'] + ' ' + df_m['Region2']
del df_m['Region1']
del df_m['Region2']

#%% Mobility Coefficients and Region Data
dict_m = df_m.set_index('Between').T.to_dict('list')
dict_r = df_r.set_index('Region').T.to_dict('dict')

#%% Testing
Tf = 180
dt = 1
r1 = RGN_SIR(dict_r['Georgia'],Tf)
r2 = RGN_SIR(dict_r['Florida'],Tf)
r3 = RGN_SIR(dict_r['Alabama'],Tf)
r4 = RGN_SIR(dict_r['California'],Tf)
r5 = RGN_SIR(dict_r['Iowa'],Tf)
r6 = RGN_SIR(dict_r['Texas'],Tf)
r7 = RGN_SIR(dict_r['New York'],Tf)
solver = SIR(r1,r2,r3,r4,r5,r6,r7)
for i in range(7):
    solver.plot_Heun(dict_m, i, dt, Tf, save = True)
