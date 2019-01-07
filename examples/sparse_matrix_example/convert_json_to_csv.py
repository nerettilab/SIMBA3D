'''
TEMPORARY SCRIPT TO READ CSV FORMAT
'''
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json
delimiter=',';
file_name='cell1_mESC_M19_a.json'
data={}
with open(file_name, 'r') as jsonfile:
    data=json.load(jsonfile)
file_name='cell1_mESC_M19_a.csv'
n=len(data['count']) # the number of entries
with open(file_name, 'w') as jsonfile:
    # write the header
    jsonfile.write('# First line:row dimension'+delimiter+'column dimensions'+delimiter+'number of entries\n') 
    jsonfile.write('# Remaining entries:row index'+delimiter+'column index'+delimiter+'interaction value\n') 
    jsonfile.write(str(data['row_dimension'])+delimiter+str(data['column_dimension'])+delimiter+str(n)+'\n')
    for ii in range(n):
       jsonfile.write(str(data['row_index'][ii])+delimiter+str(data['column_index'][ii])+delimiter+str(data['count'][ii])+'\n')
'''
#%% Use JSON to write the tasklist to file
taskfilepath=''
with open(taskfilepath, 'w') as tasklist:
    json.dump(taskss, tasklist,indent=4,sort_keys=True) # pretty print the JSON to make it human readible
    #json.dump(tasks, tasklist) # if you don't care about pretty printing the JSON just dump it this way
'''
