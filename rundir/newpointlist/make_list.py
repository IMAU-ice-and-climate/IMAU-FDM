#import modules 
import os

# Change the current working directory
os.chdir('/perm/nld3562/MSc_thesis/NPNF_Run_ERA5_george_2km/')

# Directory containing the data files
dir3 = "/ec/res4/scratch/nld3562/MSc_thesis/2km_data/era5_george/output/"

# Prefix for the data files
pre_fname1 = "ECMWF_ANT2_FDMv12AD_RACMO2.3p2_ERA5_1979-2021"

# Open old pointlist
with open('pointlist_5.txt', 'r') as file:
    old_pointlist = [int(line.strip()) for line in file]

# Create an empty list to store the numbers of missing files
listt = []

# Loop through each point from 1 to numpoints-1
for n in old_pointlist:
    numb = n
    fname = dir3 + pre_fname1 + "_1D_" + str(numb) + ".nc"
    
    # Check if the file exists
    if os.path.exists(fname):
        print(numb)  # File found, print the number of the point
    else:
        listt.append(numb)  # File not found, add the point number to the list of missing files

# Write the list of missing file numbers to a text file
with open('pointlist.txt', 'w') as filehandle:
    for listitem in listt:
        filehandle.write('%s\n' % listitem)
