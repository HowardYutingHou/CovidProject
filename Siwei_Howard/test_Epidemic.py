# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 17:20:26 2021
@author: vivid
"""

from Region import RegionEpi
from Regions import regions
import pandas as pd

import time

class Timer:
    # We create a list to store every time period Timer recorded.
    results = []
    # We initialize the attributes t1 and t2 to be 0.
    # t1 records the time when Timer starts. t2 records the time when Timer ends.
    t1 = 0
    t2 = 0
    # We initialize the attribute duration to be 0.
    duration = 0
    # We default the unit to be 1. It could also be 1/60 if the unit is minutes, or 1/3600 if it's hours.
    unit = 1
    # We make the default unit_name to be 'seconds'.
    unit_name = 'seconds'

    # The constructor method:
    def __init__(self):
        pass

    # The start method.
    def start(self):
        # If t1 is 0, we start the Timer and record the current time in t1.
        if self.t1 == 0:
            self.t1 = time.time()
        # If t1 is not 0, the Timer is already started. We return error message.
        else:
            print('Error: the Timer has already started!')

    # The end method.
    def end(self):
        # If t1 is not 0, the Timer has already started. We end the Timer, and store the current time in t2.
        if self.t1 != 0:
            self.t2 = time.time()
            # We store the time duration between t1 and t2 in duration.
            self.duration = (self.t2 - self.t1) * self.unit
            # We print out the message of the time duration.
            #print('The time taken is ' + str(self.duration) + ' ' + self.unit_name + '.')
            # We append the time duration to the list.
            self.results.append(str(self.duration) + ' ' + self.unit_name)
            # After the Timer ends, we set t1, t2, and duration back to zero.
            self.t1 = 0
            self.t2 = 0
            self.duration = 0
        # If t1 is 0, the Timer is not currently running. We return an error message.
        else:
            print('Error: the Timer is not currently running!')

    # The method to configure the Timer to display either seconds, minutes, or hours.
    def set_timer(self, input):
        # If we take 's' as parameter, we set the unit of the Timer to be seconds.
        if input == 's':
            self.unit = 1
            self.unit_name = 'seconds'
        # If we take 'm' as parameter, we set the unit of the Timer to be minutes.
        if input == 'm':
            self.unit == 1/60
            self.unit_name = 'minutes'
        # If we take 'h' as parameter, we set the unit of the Timer to be hours.
        if input == 'h':
            self.unit = 1/3600
            self.unit_name = 'hours'

    # The method to retrieve the last timer result.
    def last_result(self):
        # We simply print out the last element of the list.
        return self.results[len(self.results) - 1]



# because installation always false 
timer = Timer()

#timer.start()
# use your own files
df = pd.read_csv(r"C:\Users\Administrator\Desktop\Models\data\fake_data.csv")
df_dict = pd.read_csv(r"C:\Users\Administrator\Desktop\Models\data\fake_inter.csv")

#df = pd.read_csv(r"C:\Users\Administrator\Desktop\3_data.csv")
#df_dict = pd.read_csv(r"C:\Users\Administrator\Desktop\3_inter.csv")

# Turn this into a name dictionary
df_dict['Pair'] = df_dict['State1'] + " " + df_dict['State2']

df_dict['List'] = df_dict['Pair']
for i in range(0, df_dict.shape[0]):
    empty = []

    empty.append(df_dict['s_i'][i])
    empty.append(df_dict['s_o'][i])
    empty.append(df_dict['i_i'][i])
    empty.append(df_dict['i_o'][i])
    empty.append(df_dict['r_i'][i])
    empty.append(df_dict['r_o'][i])
    '''
    empty.append(0)
    empty.append(0)
    empty.append(0)
    empty.append(0)
    empty.append(0)
    empty.append(0)
    '''
    df_dict['List'].iloc[i] = empty

help_dict = df_dict.set_index('Pair')['List'].to_dict()

timer.start()

list_of_regions = []
for i in range(0, df.shape[0]):
    r = RegionEpi(df['S'][i], df['I'][i], df['R'][i], df['N'][i], df['beta'][i], df['gamma'][i], df['name'][i], [0,0])
    list_of_regions.append(r)

for j in list_of_regions:
    args = []
    for k in list_of_regions:
        args.append(k)
    j.set_borders(*args)

for main in range(0, len(list_of_regions) - 1):
    for con in range(main + 1, len(list_of_regions)):
        str_key = list_of_regions[main].name + " " + list_of_regions[con].name
        data = help_dict[str_key]
        list_of_regions[main].set_InterCoeffs(list_of_regions[con], data[0], data[1], data[2], data[3], data[4],
                                              data[5])

#timer.end()
#time1 = timer.last_result()

#timer.start()

input_list = []
input_list.append(100)
input_list.append(1)

for r in list_of_regions:
    input_list.append(r)

area = regions(*input_list)

print("\nThe big matrix at t=0:")
print(area.onebig_matrix)

area.Heun_solver()
# area.bruteforce_solver()
# visualize the solutions from day 0 to day T

timer.end()

area.plot_solution()

#timer.end()

print("\nColumn matrix of the du/dt, composed of s, i, r of each region at time Tf:")
print(area.column_matrix)

print("\nThe sir_over_time matrix is meant to record the s, i, r values for the regions at each time t.")
print(len(area.sir_over_time) - 1)

#print('\nTime to implement the data and initiate the separate regions is: '+time1)
print('Time to run the solve the model using the solver is: '+timer.last_result())
