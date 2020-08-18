# -*- coding: utf-8 -*-
"""
This contains the command line interface for simba3d.

This one is for the case when you either don't want to use multiprocessing,
or if you are trying to debug an issue. Doesn't provide much information about
what actually failed. This is probably somehow my fault that the errors are not
caught.

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
from simba3d.simba3d_taskrun_SOLO import main

if __name__ == "__main__":
   main()
