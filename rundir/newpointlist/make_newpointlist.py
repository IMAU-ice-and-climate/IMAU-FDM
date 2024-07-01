#%% IMPORT MODULES

import numpy as np
import pandas as pd

#%% MAKE LIST OF COMPLETED POINTS

reference_path = '/perm/rumb/code/DATA/IN_ll_FGRN055_GrIS_GIC_implicit.txt'
points_path = ''

loc_list = pd.read_csv(reference_path, header=None, names=['lon', 'lat', 'nor', 'impexp', 'cost', 'rlat', 'rlon'])
loc_ini_done = pd.read_csv('points_ini_done.txt', header=None, names=['string'])
loc_1d_done = pd.read_csv('points_1D_done.txt', header=None, names=['string'])
loc_2d_done = pd.read_csv('points_2D_done.txt', header=None, names=['string'])
loc_2ddet_done = pd.read_csv('points_2Ddet_done.txt', header=None, names=['string'])

ints_ini_done = np.zeros_like(loc_ini_done)
ints_1d_done = np.zeros_like(loc_1d_done)
ints_2d_done = np.zeros_like(loc_2d_done)
ints_2ddet_done = np.zeros_like(loc_2ddet_done)

for i in range(len(loc_ini_done['string'])):
    string_ini = loc_ini_done['string'][i]
    string_ini = string_ini[len('/ec/res4/scratch/rumb/data/output/era055/newrun/ECMWF_FGRN055_imp_run_ini_'):]
    string_ini = string_ini[:-len('.nc')]
    ints_ini_done[i] = int(string_ini)
for i in range(len(loc_1d_done['string'])):
    string_1d = loc_1d_done['string'][i]
    string_1d = string_1d[len('/ec/res4/scratch/rumb/data/output/era055/newrun/ECMWF_FGRN055_imp_run_1D_'):]
    string_1d = string_1d[:-len('.nc')]
    ints_1d_done[i] = int(string_1d)
for i in range(len(loc_2d_done['string'])):
    string_2d = loc_2d_done['string'][i]
    string_2d = string_2d[len('/ec/res4/scratch/rumb/data/output/era055/newrun/ECMWF_FGRN055_imp_run_2D_'):]
    string_2d = string_2d[:-len('.nc')]
    ints_2d_done[i] = int(string_2d)
for i in range(len(loc_2ddet_done['string'])):
    string_2ddet = loc_2ddet_done['string'][i]
    string_2ddet = string_2ddet[len('/ec/res4/scratch/rumb/data/output/era055/newrun/ECMWF_FGRN055_imp_run_2Ddetail_'):]
    string_2ddet = string_2ddet[:-len('.nc')]
    ints_2ddet_done[i] = int(string_2ddet)

ints_ini_done = np.squeeze(np.sort(ints_ini_done))
ints_1d_done = np.squeeze(np.sort(ints_1d_done))
ints_2d_done = np.squeeze(np.sort(ints_2d_done))
ints_2ddet_done = np.squeeze(np.sort(ints_2ddet_done))
ints_done = np.unique(np.concatenate([ints_ini_done, ints_1d_done, ints_2d_done, ints_2ddet_done])) # ints fully or partially done

#%% MAKE LIST POINTS STILL TO DO

ints_ini_todo = []
ints_1d_todo = []
ints_2d_todo = []
ints_2ddet_todo = []
len_list = len(loc_list)
print(f'#{len_list} in total')
for j in range(len_list):
    ipoint = j+1
    if not ipoint in ints_ini_done: ints_ini_todo.append(ipoint)
    if not ipoint in ints_1d_done: ints_1d_todo.append(ipoint)
    if not ipoint in ints_2d_done: ints_2d_todo.append(ipoint)
    if not ipoint in ints_2ddet_done: ints_2ddet_todo.append(ipoint)
    if np.mod(ipoint, 2000) == 0: print(f'#{ipoint} of {len_list} points done')

ints_todo = np.array(np.concatenate([ints_ini_todo, ints_1d_todo, ints_2d_todo, ints_2ddet_todo]))
ints_todo = np.sort(np.unique(ints_todo))

print('done')

#%% MAKE LIST OF POINTS TO REMOVE

ints_remove = []
for k in ints_todo:
    if k in ints_done:
        ints_remove.append(k)

#%% WRITE REMOVELIST.TXT FILE

count = 0
len_remove = len(ints_remove)
with open('removelist.txt','w+') as output:
    for ipoint2 in ints_remove:
        output.write(f'{ipoint2}\n')
        count += 1
        if np.mod(count, 2000) == 0: print(f'#{count} of {len_remove} points done')
print('done')

#%% WRITE POINTLIST.TXT FILE

count2 = 0
len_todo = len(ints_todo)
with open('pointlist_todo.txt','w+') as output:
    for ipoint3 in ints_todo:
        output.write(f'{ipoint3}\n')
        count2 += 1
        if np.mod(count2, 2000) == 0: print(f'#{count2} of {len_todo} points done')
print('done')

#%% WRITE NF POINTLISTS FILES

count3 = 0
for ifile in np.arange(0,10):
    count4 = 0
    nf_file_name = f'pointlist_nf_{ifile+1}.txt'
    with open(nf_file_name, 'w+') as output:
        while count4 < 300:
            ipoint4 = ints_todo[count3]
            output.write(f'{ipoint4}\n')
            count3 += 1
            count4 += 1
print('done')

#%% WRITE NP POINTLISTS FILES

count5 = 0
for ifile in np.arange(0,10):
    np_file_name = f'pointlist_np_{ifile+1}.txt'
    start_ints = int(ifile/10.*len_todo)
    end_ints = int((ifile+1)/10.*len_todo)
    with open(np_file_name, 'w+') as output:
        for ipoint5 in ints_todo[start_ints:end_ints]:
            output.write(f'{ipoint5}\n')
            count5 += 1
            if np.mod(count5, 2000) == 0: print(f'#{count5} of {len_todo} points done')
print('done')

# %%
