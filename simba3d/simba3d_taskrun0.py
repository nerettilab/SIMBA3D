# -*- coding: utf-8 -*-
"""
This contains the command line interface for simba3d.

I was asked to make command line options for all the inputs, so this is what
the command line looked like before that change took place.

Created on Thu Sep 07 06:48:42 2017

@author: Michael Rosenthal
"""

# import simplejson or json
#
# This is useful because the syntax is pythonic, it's human readable, and it
# can be read by other programs in other languages.
#
# json is potentially and older version of simplejson which is updated more
# frequently than Python
#
# A good practice is to use one or the other as a fallback.
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json
import sys
import time
import simba3d.mp_manager as mp

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

#%%
def printhelp():
    print('simba3d -c <number_of_cores> -r <json_tasklist_file>')
def main(args=None):
   if args is None:
       args=sys.argv[:]
       #print args
   inputfile = ''
   cores=1;
   ii=1;
   if len(args)==1:
       printhelp()
       sys.exit()
   while ii < len(args):
       if args[ii]== '-h':
           printhelp()
           sys.exit()
       elif args[ii]== '-c':
           ii+=1
           cores=int(args[ii] )
       elif args[ii]== '-r':
           ii+=1
           inputfile=args[ii]
       ii+=1
   print('Input file is ', inputfile)
   print('number of cores', str(cores))
   with open(inputfile,'r') as tasklist:
       tasks=json.load(tasklist)
       t= time.time()
       mp.mp_handler(tasks,cores)
       elapsed=time.time()-t
       print('Total %s seconds elapsed' % (elapsed))

if __name__ == "__main__":
   main()
