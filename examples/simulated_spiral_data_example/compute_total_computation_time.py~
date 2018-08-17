# -*- coding: utf-8 -*-
"""
Compute the total computation time from a list of simba3d output report files.

@author: Michael Rosenthal
"""
import os

import simba3d.mp_manager as mp


# specify the directory where the results are stored
report_directory='results/'

# manually enter a list of files to add coputation times
report_files =[
        'simulated_data_b7b1fe84-8b80-4bb4-9af8-8f0f9e7aa3b7_50.npz',
        'simulated_data_b7b1fe84-8b80-4bb4-9af8-8f0f9e7aa3b7_100.npz',
        'simulated_data_b7b1fe84-8b80-4bb4-9af8-8f0f9e7aa3b7_200.npz',
        'simulated_data_b7b1fe84-8b80-4bb4-9af8-8f0f9e7aa3b7_400.npz'
        ]
        

total_time=0.0     # initialize the total time  
for report_file in report_files:
  # load the summary report
  outputfile=os.path.join(report_directory,report_file) # path the file
  summary=mp.load_result(outputfile)
  print(outputfile+"\n\t"+str(summary['computation_time'])+" seconds")
  # for each report_file, load it into memory
  total_time+=summary['computation_time']

print("total time: "+str(total_time)+" seconds")
print("total time: "+str(total_time/60)+" minutes")
