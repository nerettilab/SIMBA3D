'''
TEMPORARY SCRIPT TO READ CSV FORMAT
'''
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json
delimiter=',';
file_name='cell1_mESC_M19_a.csv'
data={}
with open(file_name, 'r') as csvfile:
    # first two lines are the header
    line=csvfile.readline()
    line=csvfile.readline()
    # fist nontrivial line
    line=csvfile.readline()
    for line in csvfile:
        print(line)
        
    
