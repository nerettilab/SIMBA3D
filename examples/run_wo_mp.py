# -*- coding: utf-8 -*-
"""
Run a single task without the multiprocessing module
Created on Tue Jul 24 15:17:11 2018

@author: star
"""
# A good practice is to use one or the other as a fallback.
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json

import simba3d.mp_manager.mp

# LOAD THE TASKS

for task in tasks:
    mp.run_task(task)