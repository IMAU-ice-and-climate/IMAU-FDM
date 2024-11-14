import os
import pandas as pd

scratch = "/ec/res4/scratch/nld4814"
dir_pointlist = scratch+"/FGRN055_era055/reference/"
file_pointlist = dir_pointlist + "IN_ll_FGRN055_GrIS_GIC_implicit.txt"
dir_in = scratch+"/run-1957-2023-FGRN055-era055/output/"
filebase_in = "FGRN055_era055"

filename_out = scratch+"/run-1957-2023-FGRN055-era055/post-process/missing.txt"

#get number of points that should exist given the original pointlist
numpointlist=len(pd.read_csv(file_pointlist, header = None))

#create file if it doesn't exist and append missing point
with open(filename_out,"a+") as file:
	
	for point_i in range(numpointlist):

		filename_in = dir_in + filebase_in + "_1D_" + str(point_i+1) + ".nc"

		if not os.path.isfile(filename_in):
			file.write(f'{str(point_i+1)}\n')
			print(f'{str(point_i+1)}\n')