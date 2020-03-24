# -*- coding: utf-8 -*-
"""
Converts from SIMBA3D's JSON output to PDB format
compatible with Chimera molecular visualization software

Created on Tue March 24, 2020

@author: Nicola Neretti
"""
"""
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='sum the integers (default: find the max)')

args = parser.parse_args()
print(args.accumulate(args.integers))
"""

import json

p1 = "0.20"
p2 = "10.00"


with open("output.json") as f:
	data = json.load(f)
a = data["X_evol"][-1] # last element of X_evol contains the best solution
x = a[0]
y = a[1]
z = a[2]

np = len(x)

f = open("output.pdb", "w")
for i in range(0,np):
	s = "ATOM  %4s   CA MET A"+str(i+1)+"     "+str(i+1)+str(x[i])+" "+str(y[i])+"\t"+str(y[i])+"\t"+p1+"\t"+p2+"\n"
	f.write(s)
for i in range(0,np-1):
	s = "CONNECT\t" + str(i+1) + "\t" + str(i+2) + "\n"
	f.write(s)
f.close()

