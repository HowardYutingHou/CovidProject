# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 17:20:26 2021

@author: vivid
"""


from Region import RegionEpi
from Regions import regions
import numpy as np
import pandas as pd

# use your own files 
df = pd.read_csv(r"C:\Users\vivid\Desktop\fake_data.csv")
df_dict = pd.read_csv(r"C:\Users\vivid\Desktop\fake_inter.csv")

# Turn this into a name dictionary 
df_dict['Pair'] = df_dict['State1'] + " " + df_dict['State2']

df_dict['List'] = df_dict['Pair']
for i in range(0,df_dict.shape[0]):
    empty = []
    empty.append(df_dict['s_i'][i])
    empty.append(df_dict['s_o'][i])
    empty.append(df_dict['i_i'][i])
    empty.append(df_dict['i_o'][i])
    df_dict['List'].iloc[i] = empty 

help_dict = df_dict.set_index('Pair')['List'].to_dict()

list_of_regions = []
for i in range(0,df.shape[0]):
    r = RegionEpi(df['S'][i], df['I'][i], df['R'][i], df['N'][i], df['beta'][i] ,df['gamma'][i], df['name'][i])
    list_of_regions.append(r)

for j in list_of_regions:
    args = []
    for k in list_of_regions:
        args.append(k)
    j.set_borders(*args)


for main in range(0,len(list_of_regions)-1):
    for con in range(main+1,len(list_of_regions)):
        str_key = list_of_regions[main].name + " " + list_of_regions[con].name
        data = help_dict[str_key]
        list_of_regions[main].set_InterCoeffs(list_of_regions[con],data[0],data[1],data[2],data[3])

input_list = []
input_list.append(100)
input_list.append(1)
for r in list_of_regions:
    input_list.append(r)
area = regions(*input_list)

area.Heun_solver()
# visualize the solutions from day 0 to day T
area.plot_solution()
    
print("\nColumn matrix of the du/dt, composed of s, i, r of each region at time Tf:")
print(area.column_matrix)

print("\nThe big matrix at last iteration:")
print(area.onebig_matrix)

print("\nThe sir_over_time matrix is meant to record the s, i, r values for the regions at each time t.")
print(len(area.sir_over_time) - 1)

