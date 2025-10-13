import pandas as pd
import string
import glob
from datetime import datetime, timedelta


def create_datetime(start_date, end_date, ndays_timestep,resample_t=None):

    ##################

    # INPUTS
        # start_date: start date of dataset
        # end_date: end date of dataset
        # ndays_timestep: output frequency in days (usually 1, 10, or 30 days, specified in start_model_ccab)

    # OUTPUTS
        # date_list: array of datetime time stamps from start to end date with delta ndays_timestep

    ####################


    # Initialize an empty list
    date_list = []
    date = start_date
    
    # Loop through the range of dates and append to the list
    while date <= end_date:
        date_list.append(date)
        date += timedelta(days=ndays_timestep)

    #resample to specified frequency if desired
    if resample_t is not None:
        ts = pd.to_datetime(date_list)
        pds = pd.Series(index=ts)
        date_list = pds.resample(str(resample_t)).sum().index.to_pydatetime()

    return date_list

def set_fdm_df_time_coordinates(df, ID_num, rlat_i, rlon_i, ndays_timestep, start_date=datetime(1939,9,1), end_date=datetime(2023,12,31)):

    """

    INPUTS
        df: fdm data output (1d, 2d, 2ddetail)
        ID_num: FDM ID number (index in reference pointlist file)
        rlat_i: xx/location in grid (column 5 from reference pointlist file)
        rlon_i: yy/location in grid (column 6 from reference pointlist file)
        ndays_timestep: output frequency in days (usually 1, 10, or 30 days, specified in start_model_ccab)    
        start_date: start date of dataset
        end_date: end date of dataset

    OUTPUTS
        df: data frame with datetime time stamps as coordinate & dimension

    """
    
    date_list = create_datetime(start_date, end_date, ndays_timestep)

    df = df.isel(ind_t=slice(0,len(date_list))) # crops extraneous nans out of dataframe

    df["FDM_ID"] = ID_num # log point number

    # # you could also have "point" as a coordinate/dimension, but it refers to the same point as rlat/rlon
    # # e.g.,
    # # df = df.assign_coords({"FDM_ID":ID_num})
    # # df = df.expand_dims(dim={'FDM_ID':1})

    df = df.assign_coords(time=("ind_t",date_list)) #creates time as the coordinatedf = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t
    df = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t

    return df


def read_model_settings(dir_project, specific_point=False):

    """
    Reads first model settings file in list from an output project folder
    where variables are defined as:
    
        value ! variable_name; description

    Inputs:
        dir_project = filepath to project directory (usually in scratch)
        specific_point = False, or integer indicating which point number to to load model settings for

    Returns a dictionary with variable names as keys and values as appropriate types.
    
    # Usage
    model_settings_settings=read_model_settings(file_path_model_settings)

    for var_name, value in model_settings.items():
        print(f"{var_name}: {value}")
    """
    
    model_settings = {}


    if isinstance(specific_point, int):
        
        file_list = [f for f in glob.glob(dir_project + "ms_files/" + "*_"+str(specific_point)+".txt")]
            
        if not file_list:
            print("Point " + str(specific_point) + " not found in " + dir_project + "ms_files/. Using the first model settings file in the directory.")
            specific_point = False
    
    if not specific_point:
        file_list = [f for f in glob.glob(dir_project + "ms_files/" + "*.txt")]
        
    else:
        print("Specific point option, " + str(specific_point) + " not recognized.")
        return None

    f = file_list[0]
    n = 0 # to deal with duplicate names
    
    with open(f, 'r') as file:
        for line in file:
            line = line.strip()
            
            # Skip empty lines and comment-only lines
            if not line or line.startswith('!'):
                continue
            
            # Check if line contains a value and variable definition
            if '!' in line:
                # Split on the first exclamation mark
                parts = line.split('!', 1)
                value_str = parts[0].strip()
                comment_part = parts[1].strip()
                
                # Extract variable name (first word after !)
                if comment_part:
                    var_name = comment_part.split()[0].rstrip(string.punctuation)
                    if n==0:
                        var_name = "nyears_total"
                        n=n+1
                    elif n==1:
                        var_name = "nyears_spinup"
                        n=n+1
                    else:
                        n=n+1
                    
                    # Try to convert value to appropriate type
                    try:
                        # Try integer first
                        if '.' not in value_str:
                            value = int(value_str)
                        else:
                            # Try float
                            value = float(value_str)
                    except ValueError:
                        # Keep as string if conversion fails
                        value = value_str
                    
                    model_settings[var_name] = value
    
    return model_settings