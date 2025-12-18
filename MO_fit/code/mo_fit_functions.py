import pandas as pd
import string
import glob
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from process_sumup_cores import def_models
from scipy.optimize import curve_fit
from scipy.optimize import fsolve, brentq

def prep_2D_dataset(ds):
    """
    Prepares a 2D FDM dataset by removing all time steps where all variables are NaN across all layers.
    
    Parameters:
    ds (xarray.Dataset): The input dataset with dimensions 'ind_t' (time) and 'layer' (depth).
    
    Returns:
    xarray.Dataset: The cropped dataset with only valid time steps.
    """

    # Method 1: If you have a Dataset with multiple variables
    if isinstance(ds, xr.Dataset):
        # Check where ALL variables are NaN across all layers for each ind_t
        # Stack all data variables together and check for all-NaN time steps
        all_nan_mask = np.all([
            ds[var].isnull().all(dim='layer') 
            for var in ds.data_vars
        ], axis=0)
        
        # Find the first index where all variables are NaN
        valid_indices = np.where(~all_nan_mask)[0]
        
        if len(valid_indices) > 0:
            last_valid_idx = valid_indices[-1] + 1
            ds_cropped = ds.isel(ind_t=slice(0, last_valid_idx))
        else:
            # Handle case where all data is NaN
            ds_cropped = ds.isel(ind_t=slice(0, 0))

    return ds_cropped

def load_2d_file(project_dir, file_path, clean_2D=False):

    path_len = len(f"{project_dir}/output/FGRN055_era055_2D")

    df = xr.open_dataset(file_path)

    if clean_2D:
        df = prep_2D_dataset(df)

    df["point_index"] = int(file_path[path_len:-3]) # add fdm index to imported dataset

    ndays_timestep = 30

    # add timestamp as coordinate & dimension
    date_list = create_datetime(datetime(1939,9,1), datetime(2023,12,31), ndays_timestep, resample_t=None)
    date_list = date_list[1:]

    #replace indexed coordinates with timeseries

    df = df.assign_coords(time=("ind_t",date_list)) #creates time as the coordinatedf = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t
    df = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t

    return df

def load_1d_file(project_dir, file_name, clean_1D=False):

    path_len = len(f"{project_dir}/output/FGRN055_era055_1D")

    df = xr.open_dataset(file_name)

    df["point_index"] = int(file_name[path_len:-3]) # add fdm index to imported dataset

    ndays_timestep = 1
    
    # add timestamp as coordinate & dimension
    date_list = create_datetime(datetime(1939,9,1), datetime(2023,12,31), ndays_timestep, resample_t=None)
    

    #replace indexed coordinates with timeseries
    try:
        df = df.assign_coords(time=("ind_t",date_list)) #creates time as the coordinatedf = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t
        df = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t
    except:
        date_list = date_list[1:]
        df = df.assign_coords(time=("ind_t",date_list)) #creates time as the coordinatedf = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t
        df = df.swap_dims({"ind_t":"time"}) #swaps time as the dimension from ind_t

    return df

def match_core_to_fdm(data_df, core_match, avg_yearly=True):

    """
    Matches a core to the corresponding FDM data based on year and returns the selected data.
    Parameters:
    data_df (xarray.Dataset): The FDM dataset with time dimension.
    core_match (pandas.Series): The core information including the year to match.
    avg_yearly (bool): If True, averages data over the entire year; if False, selects data at the end of the year.
    """

    # swap ind_t for real(ish) timesteps
    if avg_yearly:
        data_df_sel = data_df.sel(time=slice(datetime(core_match.year,1,1), datetime(core_match.year,12,31))).mean(dim='time')
    else:
        data_df_sel = data_df.sel(time=datetime(core_match.year,12,31))
    
    data_df_sel = data_df_sel.swap_dims({"layer":"depth"})

    # drops all depths where depth is nan
    data_df_sel = data_df_sel.dropna(dim='depth', how='all')

    # reverses so that the shallowest values are early in the array and the deeper values are later in the array; helps with curve fitting and readability
    data_df_sel = data_df_sel.isel({dim: slice(None, None, -1) for dim in data_df_sel.dims})

    return data_df_sel

def match_fdm_to_year(data_df, year, avg_yearly=True):

    """
    Matches a core to the corresponding FDM data based on year and returns the selected data.
    Parameters:
    data_df (xarray.Dataset): The FDM dataset with time dimension.
    core_match (pandas.Series): The core information including the year to match.
    avg_yearly (bool): If True, averages data over the entire year; if False, selects data at the end of the year.
    """

    # swap ind_t for real(ish) timesteps
    if avg_yearly:
        data_df_sel = data_df.sel(time=slice(datetime(year,1,1), datetime(year,12,31))).mean(dim='time')
    else:
        data_df_sel = data_df.sel(time=datetime(year,12,31))
    
    data_df_sel = data_df_sel.swap_dims({"layer":"depth"})

    # drops all depths where depth is nan
    data_df_sel = data_df_sel.dropna(dim='depth', how='all')

    # reverses so that the shallowest values are early in the array and the deeper values are later in the array; helps with curve fitting and readability
    data_df_sel = data_df_sel.isel({dim: slice(None, None, -1) for dim in data_df_sel.dims})

    return data_df_sel

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

def prep_interim_dataset(df_fdmv1p2, var):

    # find true length of data and crop extra values at end are nans
    data_len = len(df_fdmv1p2[var].isel(rlat=300,rlon=150).values[df_fdmv1p2[var].isel(rlat=100,rlon=100).notnull().values])
    df_fdmv1p2 = df_fdmv1p2.isel(time=slice(0,data_len)) 
    
    start_date_1p2 = datetime(1957,10,1)
    end_date_1p2 = datetime(2020,12,31) #in days
    ndays_timestep = 30
    
    #could be off by 1 day
    #date_list = create_datetime(start_date_1p2, end_date_1p2, ndays_timestep,resample_t=None)
    date_list_resampled_1p2 = create_datetime(start_date_1p2, end_date_1p2, ndays_timestep)
    date_list_resampled_1p2 = date_list_resampled_1p2[1:]
    
    df_fdmv1p2["time"] = date_list_resampled_1p2

    file_path_pointlist = "../../reference/FGRN055/IN_ll_FGRN055.txt"
    
    pointlist_df = pd.read_csv(file_path_pointlist,names=["longitude","latitude","rlat","rlon"],usecols=[0,1,5,6])
    
    df_fdmv1p2["FDM_ID"] = (["rlat","rlon"],np.full([len(df_fdmv1p2.rlat),len(df_fdmv1p2.rlon)],np.nan))

    for fdm_i in range(0,len(pointlist_df)):

        try:
            rlat_i = pointlist_df["rlat"].iloc[fdm_i].astype(int)		  	# Obtain x location		
            rlon_i = pointlist_df["rlon"].iloc[fdm_i].astype(int)		  	# Obtain y location 
            
            df_fdmv1p2["FDM_ID"].loc[dict(rlat=rlat_i, rlon=rlon_i)] = int(fdm_i+1)
        except:
            print(str(fdm_i))

    df_fdmv1p2 = df_fdmv1p2.set_coords("FDM_ID")
    
    if "zs" in list(df_fdmv1p2.keys()):
        df_fdmv1p2 = df_fdmv1p2.rename({"zs":"h_surf"})

    return df_fdmv1p2

def find_FDM_depths(merged_df, file_list, project_dir, target_densities, clean_2D=False, do_plot=False, avg_yearly=True):

    all_depths = {}

    models = def_models()

    for file_i, file_name in enumerate(file_list):
        
        # load data and find core match
        try:

            data_df = load_2d_file(project_dir, file_name, clean_2D)

        # file doesn't exist or is corrupt (former shouldn't happen, latter does occcasionally - need to rerun FDM point if so)    
        except Exception as e:
            print(f"Failed to load file {file_i} - {file_name}: {e}")
            continue

        # find core(s) closest to this fdm point
        core_matches = merged_df[merged_df["FDM_point_index"] == data_df.point_index.values]

        if len(core_matches) == 0:
            print(f"No core match for FDM point {data_df.point_index.values} in file {file_i} - {file_name}")
            continue

        # often more than one core is near the same FDM point
        for core_i in range(len(core_matches)):
            
            core_match = core_matches.iloc[core_i]

            data_df_sel = match_core_to_fdm(data_df, core_match, avg_yearly) # selects the correct date for each core
            
            density = data_df_sel.dens.values
            depth = data_df_sel.depth.values

            depth_smooth = np.linspace(depth.min(), depth.max(), 200)

            fitted_params = {}

            best_r2 = -np.inf
            best_model = None

            for idx, (func, name, expression, color, func_min) in enumerate(models):

                str_point = str(data_df_sel["point_index"].values)

                try:
                    str_core = core_match["core_name"]
                    title_core = f'FDM: {str_point} Core: {str_core}'
                except:
                    str_core = str(core_match.name)
                    title_core = f'FDM: {str_point} Core: {str_core}'

                try:
                    # Fit the model
                    popt, pcov = curve_fit(func, depth, density)
                    
                    # Store parameters with a simple key
                    model_name = name.split(":")[0].strip().lower().replace(" ", "_")
                    fitted_params[model_name] = popt
                    
                    # Calculate fitted values
                    density_fit = func(depth, *popt)
                    density_smooth = func(depth_smooth, *popt)
                    
                    # Calculate R²
                    residuals = density - density_fit
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((density - np.mean(density))**2)
                    r_squared = 1 - (ss_res / ss_tot)

                    if r_squared > best_r2:
                        best_r2 = r_squared
                        best_model = (func, popt, name, expression, round(best_r2,3), func_min)
                    
                    # if do_plot:

                        # fig, axes = plt.subplots(figsize=(5, 5))
                            
                        # # Plot data and fit (top row)
                        # ax1 = axes
                        # ax1.scatter(depth, density, alpha=0.6, label='Model', s=30)
                        # ax1.plot(depth_smooth, density_smooth, color=color, linewidth=2, label='Fit')
                        # ax1.plot(core_match['depth_to_550'], 550, 'm*', label='Core 550 kg/m³', markersize=8)
                        # ax1.plot(core_match['depth_to_830'], 830, 'm*', label='Core 830 kg/m³', markersize=8)
                        # ax1.set_xlabel('Depth')
                        # ax1.set_ylabel('Density')
                        # ax1.legend()
                        # ax1.grid(True, alpha=0.3)
                        # ax1.set_title(f'{title_core}\n{name.split(":")[0]}\nR² = {r_squared:.4f}')
                        
                    #     # Plot residuals (bottom row)
                    #     ax2 = axes[1, idx]
                    #     ax2.scatter(depth, residuals, alpha=0.6, color=color, s=30)
                    #     ax2.axhline(y=0, color='black', linestyle='--', linewidth=1)
                    #     ax2.set_xlabel('Depth')
                    #     ax2.set_ylabel('Residuals')
                    #     ax2.grid(True, alpha=0.3)
                    #     ax2.set_title('Residuals')
                    
                except Exception as e:
                    print(f"\n{name}: Failed to fit - {e}")
                    
            #         if do_plot:
                        
            #             fig, ax1 = plt.subplots(1, figsize=(5,5))

            #             ax1.scatter(depth, density, alpha=0.6, label='Model', s=30)
            #             ax1.plot(core_match['depth_to_550'], 550, 'm*', label='Core 550 kg/m³', markersize=8)
            #             ax1.plot(core_match['depth_to_830'], 830, 'm*', label='Core 830 kg/m³', markersize=8)
            #             ax1.set_xlabel('Depth')
            #             ax1.set_ylabel('Density')
            #             ax1.legend()
            #             ax1.grid(True, alpha=0.3)
            #             ax1.set_title(f'{title_core}')
                    
            # if do_plot:
            #     plt.tight_layout()
            #     plt.show()

            target_depths = {}
            failed = 0

            for rho in target_densities:
                #print(f"\nFinding depth for target density: {rho} kg/m³ using a {best_model[2]} fit")

                if best_model[4] < 0.9:
                    #print(f"R² too low ({best_model[4]}) to reliably estimate depth for density {rho} kg/m³")
                    failed = int(1)
                    depth_at_rho = np.nan
                
                else:

                    try: 
                        depth_at_rho = brentq(
                            best_model[5],
                            depth.min(),
                            depth.max(),
                            args=(*best_model[1], rho)
                        )


                    except ValueError as e:

                        depth_at_rho = np.nan
                        failed = int(1)

                target_depths[rho] = {"depth":round(depth_at_rho,2), "r2":float(best_model[4]), "model_name":best_model[2], "failed":failed}

            if (do_plot) and (failed==0):
                    
                fig, ax1 = plt.subplots(1, figsize=(5,5))

                ax1.scatter(depth, density, alpha=0.6, color='grey', label='Model', s=30)
                ax1.plot(core_match['depth_to_550'], 550, 'm*', label='Core 550 kg/m³', markersize=8)
                ax1.plot(core_match['depth_to_830'], 830, 'm*', label='Core 830 kg/m³', markersize=8)
                ax1.plot(target_depths[550]["depth"], 550, 'g^', label='Fit 550 kg/m³', markersize=8)
                ax1.plot(target_depths[830]["depth"], 830, 'g^', label='Fit 830 kg/m³', markersize=8)
                ax1.set_xlabel('Depth')
                ax1.set_ylabel('Density')
                ax1.legend()
                ax1.grid(True, alpha=0.3)
                ax1.set_title(f'{title_core}')
                
                plt.tight_layout()

                plt.save_figure("../figures/")

            all_depths[int(core_match.name)] = target_depths


    # add depths to merged df
    for i, row_i in enumerate(all_depths):
        merged_df.loc[row_i, "FDM_depth_to_550"] = all_depths[row_i][550]["depth"]
        merged_df.loc[row_i, "FDM_depth_to_830"] = all_depths[row_i][830]["depth"]
        merged_df.loc[row_i, "Fit_R2"] = all_depths[row_i][550]["r2"]
        merged_df.loc[row_i, "Fit_model"] = all_depths[row_i][550]["model_name"]
        merged_df.loc[row_i, "Fit_failed"] = all_depths[row_i][550]["failed"]

    #removes all failed or nan fits
    merged_df = merged_df[merged_df["Fit_failed"]==0]
    merged_df = merged_df[~merged_df["Fit_failed"].isna()]

    return merged_df

def add_vars_to_merged_df(merged_df, forcing_dir="/home/nld4814/scratch/FGRN055_era055/input/averages/"):
    
    mask = xr.open_dataset("../../reference/FGRN055/FGRN055_Masks.nc")

    for var in ["precip", "ff10m", "tskin", "snowfall"]:

        df = xr.open_dataset(f"{forcing_dir}{var}_FGRN055_era055-1939_1940-1970_ave.nc")
        df = df.squeeze()
        df = df.drop_vars(['assigned','block1', 'block2', 'height','time'])
        df["rlat"] = mask["rlat"]
        df["rlon"] = mask["rlon"]

        for idx in merged_df.index: 

            rlat_i =int(merged_df.loc[idx,"FDM_rlat"])
            rlon_i = int(merged_df.loc[idx, "FDM_rlon"])

            df_point = df[f"{var}"].sel(rlat=rlat_i, rlon=rlon_i)
            
            merged_df.loc[idx,f"{var}"] = float(df_point.values)

        df.close()

    return merged_df

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

def set_fdm_df_time_coordinates(df, ID_num, ndays_timestep, start_date=datetime(1939,9,1), end_date=datetime(2023,12,31)):

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


    if not isinstance(specific_point, bool):
        
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

def match_fdm_points_to_cores(core_locations, fdm_pointlist_df, save_pointlist=False, output_dir="../data/merged/", do_plot=False):

    xs = []; distances = []  # distance between the pair of points
    
    for point in core_locations:

        assert len(point) == 2, "``points`` should be a tuple or list of tuples (lat, lon)"

        p_lat, p_lon = point
        # Find absolute difference between requested point and the grid coordinates.
        abslat = np.abs(fdm_pointlist_df.latitude - p_lat)
        abslon = np.abs(fdm_pointlist_df.longitude - p_lon)

        # Create grid of the maximum values of the two absolute grids
        c = np.maximum(abslon, abslat)

        # Find location where lat/lon minimum absolute value intersects
        x = np.where(c == np.min(c))[0][0]
        xs.append(x)

        # Matched Grid lat/lon
        g_lat = fdm_pointlist_df.iloc[x,:].latitude
        g_lon = fdm_pointlist_df.iloc[x,:].longitude

        R = 6373.0  # approximate radius of earth in km

        lat1 = np.deg2rad(p_lat); lon1 = np.deg2rad(p_lon)
        lat2 = np.deg2rad(g_lat); lon2 = np.deg2rad(g_lon)
        dlon = lon2 - lon1; dlat = lat2 - lat1

        a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

        distance = R * c
        distances.append(distance)

    #pointlist_df.loc[xs].plot.scatter(x='longitude', y='latitude', title='FDM points nearest to cores', xlabel='Longitude', ylabel='Latitude',figsize=(4,5))

    if save_pointlist:
        output_path = f"{output_dir}pointlist_from_cores.csv"
        pointlist_near_cores = fdm_pointlist_df.loc[xs].index.values
        np.savetxt(output_path, pointlist_near_cores, fmt='%d')

    if do_plot:
        
        fig, ax = plt.subplots()

        ax.scatter(core_locations[:,1], core_locations[:,0], label='Cores')
        fdm_pointlist_df.loc[xs].plot.scatter(x='longitude', y='latitude', color='red', marker='x', alpha=0.7, title='FDM points nearest to cores', xlabel='Longitude', ylabel='Latitude',ax=ax, label='FDM points')
        plt.legend()

    return xs


#duplicate with def_models
def def_MO_models():
    # Define fitting functions for MO
    def mo_logarithmic(acav, a, b):
        """MO = a + b*log(acav)"""
        return a + b * np.log(acav)

    def mo_powerlog(acav, a, b, c):
        """MO = a * acav^b + c"""
        return a * (acav**b) + c
    
    models = [
        mo_powerlog,
        mo_logarithmic
    ]    

    return models



def do_MO_fit(df_clean, var="snowfall"):

    """Performs MO fitting for both 550 and 830 kg/m³ density regimes.
    
    Parameters:
    df_clean (pandas.DataFrame): DataFrame containing core and modeled depths to 550 and 830.
    df_type (str): which model run to use ("interim", "era055", "MO").
    Returns:
    Dictionary of fit results for both density regimes.
    """

    df_clean.loc[:,'MO_optimal_550'] = df_clean[f'FDM_depth_to_550_MO']/df_clean['depth_to_550']
    df_clean.loc[:,'MO_optimal_830'] = (df_clean[f'FDM_depth_to_830_MO']-df_clean[f'FDM_depth_to_550_MO'])/(df_clean['depth_to_830']-df_clean['depth_to_550'])
    

    # Remove unrealistic MO values (e.g., < 0.1 or > 5)
    df_clean = df_clean[(df_clean['MO_optimal_550'] > 0.25) & (df_clean['MO_optimal_550'] < 1)]

    mo_powerlog, mo_logarithmic = def_MO_models()

    # Perform fits for 830 density regime
    print("\n" + "="*60)
    print("FITTING MO FOR 550 kg/m³ DENSITY REGIME")
    print("="*60)

    spy = 365.25*24*3600
    var_avg = df_clean[var]
    
    df_clean_830 = df_clean[~df_clean['MO_optimal_830'].isna()]
    var_avg_830 = df_clean_830[var]

    if (var=="precip") or (var=="snowfall"):
        var_avg = var_avg*spy
        var_avg_830 = var_avg_830*spy

    # Logarithmic fit for 550
    popt_log_550, pcov_log_550 = curve_fit(mo_logarithmic, var_avg, df_clean['MO_optimal_550'], maxfev = 2000)
    a_log_550, b_log_550 = popt_log_550
    MO_pred_log_550 = mo_logarithmic(var_avg, *popt_log_550)
    r2_log_550 = 1 - np.sum((df_clean['MO_optimal_550'] - MO_pred_log_550)**2) / np.sum((df_clean['MO_optimal_550'] - df_clean['MO_optimal_550'].mean())**2)
    rmse_log_550 = np.sqrt(np.mean((df_clean['MO_optimal_550'] - MO_pred_log_550)**2))

    print(f"\nLogarithmic fit: MO = {a_log_550:.4f} + {b_log_550:.4f}*log({var})")
    print(f"  R² = {r2_log_550:.4f}")
    print(f"  RMSE = {rmse_log_550:.4f}")
    if (var=="snowfall") or (var=="precip"):
        print(f"Current FGRN fit: MO = 0.6688 + 0.0048*log(acav)")

    # Perform fits for 830 density regime
    print("\n" + "="*60)
    print("FITTING MO FOR 830 kg/m³ DENSITY REGIME")
    print("="*60)

    # Logarithmic fit for 830
    popt_log_830, pcov_log_830 = curve_fit(mo_logarithmic, var_avg_830, df_clean_830['MO_optimal_830'][~df_clean['MO_optimal_830'].isna()], maxfev = 2000)
    a_log_830, b_log_830 = popt_log_830
    MO_pred_log_830 = mo_logarithmic(var_avg_830, *popt_log_830)
    r2_log_830 = 1 - np.sum((df_clean_830['MO_optimal_830'] - MO_pred_log_830)**2) / np.sum((df_clean_830['MO_optimal_830'] - df_clean_830['MO_optimal_830'].mean())**2)
    rmse_log_830 = np.sqrt(np.mean((df_clean_830['MO_optimal_830'] - MO_pred_log_830)**2))

    print(f"\nLogarithmic fit: MO = {a_log_830:.4f} + {b_log_830:.4f}*log({var})")
    print(f"  R² = {r2_log_830:.4f}")
    print(f"  RMSE = {rmse_log_830:.4f}")
    if (var=="snowfall") or (var=="precip"):
        print(f"Current FGRN fit: MO = 1.7465 - 0.2045*log(acav)")

    # Powerlog fit for 830
    try: 
        popt_pow_830, pcov_pow_830 = curve_fit(mo_powerlog, var_avg_830, df_clean_830['MO_optimal_830'], p0=[1, -0.6, 1], maxfev = 2000)
        a_pow_830, b_pow_830, c_pow_830 = popt_pow_830
        MO_pred_pow_830 = mo_powerlog(var_avg_830, *popt_pow_830)
        r2_pow_830 = 1 - np.sum((df_clean_830['MO_optimal_830'] - MO_pred_pow_830)**2) / np.sum((df_clean_830['MO_optimal_830'] - df_clean_830['MO_optimal_830'].mean())**2)
        rmse_pow_830 = np.sqrt(np.mean((df_clean_830['MO_optimal_830'] - MO_pred_pow_830)**2))

        print(f"\nPowerlog fit: MO = {a_pow_830:.4f} * {var}^{b_pow_830:.4f} + {c_pow_830:.4f}")
        print(f"  R² = {r2_pow_830:.4f}")
        print(f"  RMSE = {rmse_pow_830:.4f}")
        if (var=="snowfall") or (var=="precip"):
            print(f"Current ANT27 powerlog: MO = 6.387 * acav^(-0.477) + 0.195")
    
    except Exception as e:
        print(f"\nPowerlog fit: Failed to fit - {e}")
        popt_pow_830, MO_pred_pow_830, r2_pow_830, rmse_pow_830 = (None, None, None, None)

    fit_results = {
        'log_550': (popt_log_550, MO_pred_log_550, r2_log_550, rmse_log_550),
        'log_830': (popt_log_830, MO_pred_log_830, r2_log_830, rmse_log_830),
        'pow_830': (popt_pow_830, MO_pred_pow_830, r2_pow_830, rmse_pow_830)
        }
    
    return df_clean, fit_results

def plot_MO_fits(df_clean, fit_results, var="snowfall"):

    mo_powerlog, mo_logarithmic = def_MO_models()
    df_clean_830 = df_clean[~df_clean['MO_optimal_830'].isna()].copy()

    popt_log_550, MO_pred_log_550, r2_log_550, rmse_log_550 = fit_results['log_550']
    popt_log_830, MO_pred_log_830, r2_log_830, rmse_log_830= fit_results['log_830']
    popt_pow_830, MO_pred_pow_830, r2_pow_830, rmse_pow_830 = fit_results['pow_830']
     
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    spy = 3600*24*365.25
    var_avg = df_clean[var]
    var_avg_830 = df_clean_830[var]

    if (var=="precip") or (var=="snowfall"):
        var_avg = var_avg*spy
        var_avg_830 = var_avg_830*spy
        var_range = np.linspace(90,1200, 200)
        xlabel = "Accumulation (kg/m²/yr)"
    elif var=="ff10m":
        var_range = np.linspace(6,9,200)
        xlabel = "FF10m (m/s)"

    # Generate smooth curves for plotting

    # Current FGRN fits for comparison
    if (var=="precip") or (var=="snowfall"):
        MO_current_550 = 0.6688 + 0.0048 * np.log(var_range)
        MO_current_830 = 1.7465 - 0.2045 * np.log(var_range)
        MO_ant_current_830 = 6.387 * (var_range)**(-0.477) + 0.195
    else: 
        MO_current_550 = None
        MO_currrent_830 = None
        MO_ant_current_830 = None

    high_avac_df = df_clean[df_clean["precip"]>3.25e-5].copy()
    high_avac_df.loc[:,'MO_optimal_550'] = high_avac_df[f'FDM_depth_to_550_MO']/high_avac_df['depth_to_550'] 

    print("Average MO at 550 is: ", np.mean(df_clean['MO_optimal_550']))
    # Plot 550 results
    ax1 = axes[0, 0]
    ax1.scatter(var_avg, df_clean['MO_optimal_550'], alpha=0.5, s=20, label='Data', color='gray')
    ax1.plot(var_range, mo_logarithmic(var_range, *popt_log_550), 'b-', linewidth=2, 
            label=f'Log fit (R²={r2_log_550:.3f})')
    if var=="snowfall":
        ax1.plot(var_range, MO_current_550, 'k--', linewidth=1.5, alpha=0.7, label='Current FGRN')
        ax1.scatter(high_avac_df[var]*spy, high_avac_df['MO_optimal_550'], facecolors='orange', edgecolors='orange', s=50, label='High avac (>3.25e-5)')
    
    ax1.set_xlabel(xlabel, fontsize=11)
    ax1.set_ylabel('MO', fontsize=11)
    ax1.set_title('MO Fits for 550 kg/m³ Density', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')


    # Plot 830 results
    ax2 = axes[0, 1]
    ax2.scatter(var_avg_830, df_clean_830['MO_optimal_830'], alpha=0.5, s=20, label='Data', color='gray')
    ax2.plot(var_range, mo_logarithmic(var_range, *popt_log_830), 'b-', linewidth=2, 
            label=f'Log fit (R²={r2_log_830:.3f})')
    if popt_pow_830 is not None:
        ax2.plot(var_range, mo_powerlog(var_range, *popt_pow_830), 'r-', linewidth=2, 
                label=f'Power fit (R²={r2_pow_830:.3f})')
    if var=="snowfall":
        ax2.plot(var_range, MO_current_830, 'k--', linewidth=1.5, alpha=0.7, label='Current FGRN')
        ax2.plot(var_range, MO_ant_current_830, 'k:', linewidth=1.5, alpha=0.7, label='Current ANT27')
    ax2.set_xlabel(xlabel, fontsize=11)
    ax2.set_ylabel('MO', fontsize=11)
    ax2.set_title('MO Fits for 830 kg/m³ Density', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')

    # Residual plots for 550
    ax3 = axes[1, 0]
    ax3.scatter(var_avg, df_clean['MO_optimal_550'] - MO_pred_log_550, 
            alpha=0.5, s=15, label='Log residuals', color='blue')
    #ax3.scatter(df_clean[var']*spy, df_clean['MO_optimal_550'] - MO_pred_pow_550, 
    #           alpha=0.5, s=15, label='Power residuals', color='red')
    ax3.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax3.set_xlabel(xlabel, fontsize=11)
    ax3.set_ylabel('Residuals', fontsize=11)
    ax3.set_title('Residuals for 550 kg/m³', fontsize=12)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    # Residual plots for 830
    ax4 = axes[1, 1]
    ax4.scatter(var_avg_830, df_clean_830['MO_optimal_830'] - MO_pred_log_830, 
            alpha=0.5, s=15, label='Log residuals', color='blue')
    if popt_pow_830 is not None:
        ax4.scatter(var_avg_830, df_clean_830['MO_optimal_830'] - MO_pred_pow_830, 
            alpha=0.5, s=15, label='Power residuals', color='red')
    ax4.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax4.set_xlabel(xlabel, fontsize=11)
    ax4.set_ylabel('Residuals', fontsize=11)
    ax4.set_title('Residuals for 830 kg/m³', fontsize=12)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
    