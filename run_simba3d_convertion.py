# -*- coding: utf-8 -*-
"""
This contains the command line interface for simba3d.

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
from simba3d.convertion_tool import main

if __name__ == "__main__":
   main()
