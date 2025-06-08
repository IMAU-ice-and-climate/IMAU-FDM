import numpy as np
import pandas as pd

username="nld4814"
project_name="FGRN055-era055_1939-2023"
#%% MAKE LIST OF COMPLETED POINTS

submission_number="2"
reference_path = "/ec/res4/scratch/"+username+"/"+project_name+"/pointlist_"+submission_number+".txt"
restart_path = "/ec/res4/scratch/"+username+"/restart/"+project_name+"/"
output_path = "/ec/res4/scratch/"+username+"/"+project_name+"/output/"
points_path = "/perm/nld4814/code/IMAU-FDM/rundir/newpointlist/"

restart_name = "FGRN055_era055_restart_from_spinup_"
output_name = "FGRN055_era055_"

loc_list = pd.read_csv(reference_path,header=None, names=['string'])
loc_1d_done = pd.read_csv(points_path+'points_1D_done.txt', header=None, names=['string'])
loc_2d_done = pd.read_csv(points_path+'points_2D_done.txt', header=None, names=['string'])
loc_2ddet_done = pd.read_csv(points_path+'points_2Ddet_done.txt', header=None, names=['string'])

ints_1d_done = []
ints_2d_done = []
ints_2ddet_done = []

for i in range(len(loc_1d_done['string'])):
    string_1d = loc_1d_done['string'][i]
    string_1d = string_1d[len(output_path+output_name+"1D_"):]
    string_1d = string_1d[:-len('.nc')]
    if string_1d:
        ints_1d_done.append(int(string_1d))
for i in range(len(loc_2d_done['string'])):
    string_2d = loc_2d_done['string'][i]
    string_2d = string_2d[len(output_path+output_name+"2D_"):]
    string_2d = string_2d[:-len('.nc')]
    if string_2d:
        ints_2d_done.append(int(string_2d))
for i in range(len(loc_2ddet_done['string'])):
    string_2ddet = loc_2ddet_done['string'][i]
    string_2ddet = string_2ddet[len(output_path+output_name+"2Ddetail_"):]
    string_2ddet = string_2ddet[:-len('.nc')]
    if string_2ddet:
        ints_2ddet_done.append(int(string_2ddet))

ints_1d_done = np.sort(ints_1d_done)
ints_2d_done = np.sort(ints_2d_done)
ints_2ddet_done = np.sort(ints_2ddet_done)
ints_done = np.unique(np.concatenate([ints_1d_done, ints_2d_done, ints_2ddet_done])) # ints fully or partially done

#%% MAKE LIST POINTS STILL TO RUN

ints_1d_todo = []
ints_2d_todo = []
ints_2ddet_todo = []
ints_todo=[]
len_list = len(loc_list)

print(f'#{len_list} in total')
for j in range(len_list):
    ipoint = j+1
    if not ipoint in ints_1d_done: ints_1d_todo.append(ipoint)
    if not ipoint in ints_2d_done: ints_2d_todo.append(ipoint)
    if not ipoint in ints_2ddet_done: ints_2ddet_todo.append(ipoint)
    if not ipoint in ints_done: ints_todo.append(ipoint)
    if np.mod(ipoint, 2000) == 0: print(f'#{ipoint} of {len_list} points done')

print('done')

#%% WRITE POINTLIST.TXT FILE

count2 = 0
len_todo = len(ints_todo)
with open(project_path+"pointlist_todo.txt",'w+') as output:
    for ipoint3 in ints_todo:
        
        output.write(f'{ipoint3}\n')
        count2 += 1
        if np.mod(count2, 2000) == 0: print(f'#{count2} of {len_todo} points done')
print('done')
