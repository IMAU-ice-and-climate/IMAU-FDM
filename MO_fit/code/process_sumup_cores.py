import pandas as pd
import numpy as np
import os
import warnings

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.optimize import brentq

warnings.filterwarnings('ignore')


"""
Extracting 550 and 830 kg/m3 depths from the SUMup Database

This script extracts core data from the SUMup database and finds depths 
to specific density thresholds (550 and 830 kg/m³) through interpolation.

The 2024 SUMup database can be downloaded from:
"https://arcticdata.io/catalog/view/doi:10.18739/A2M61BR5M"

Extracted from SUMup_2024_density_csv.zip:
- SUMup_2024_density_greenland.csv
- SUMup_2024_density_antarctica.csv  
- SUMup_2024_density_profile_names.tsv
- SUMup_2024_density_references.tsv

Authors: Claude AI Assistant and Elizabeth Case
Date: September 2025
"""

def download_sumup_data_instructions():
    """Print instructions for downloading SUMup data"""
    print("""
    To use this script, you need to download the SUMup 2023 dataset:
    
    1. Go to: https://arcticdata.io/catalog/view/doi:10.18739/A2M61BR5M
    2. Download the following CSV files:
       - SUMup_2024_density_csv.zip
            - SUMup_2024_density_greenland.csv
            - SUMup_2024_density_antarctica.csv
            - SUMup_2024_density_profile_names.tsv
            - SUMup_2024_density_references.tsv
    
    3. Place them in a folder and update the 'data_path' variable below.
    """)

def load_sumup_data(data_path,region='Greenland'):
    """
    Load SUMup density data and metadata
    
    Parameters:
    data_path (str): Path to folder containing SUMup CSV files
    region: 'Greenland' or 'Antarctica' ('TODO: add 'both' option')    
    Returns:
    tuple: (combined_density_df, profile_names_df, references_df)
    """
    
    print("Loading SUMup data files...")

    match region.lower():
        case 'greenland':
            print("Loading Greenland data only")
            # Load density data
            file = os.path.join(data_path, 'SUMup_2024_density_greenland.csv')
        case 'antarctica':
            print("Loading Antarctica data only")
            file = os.path.join(data_path, 'SUMup_2024_density_antarctica.csv')
        # case 'both':
        #     print("Loading both Greenland and Antarctica data")
        #     greenland_file = os.path.join(data_path, 'SUMup_2024_density_greenland.csv')
        #     antarctica_file = os.path.join(data_path, 'SUMup_2024_density_antarctica.csv')
        #     files = [greenland_file, antarctica_file]
        case _:
            raise ValueError("Invalid region specified. Choose 'Greenland', 'Antarctica'")#, or 'both'.")
    
    if not os.path.exists(file):
        raise FileNotFoundError(file)
        
    # Read density data
    df = pd.read_csv(file)
    df['region'] = region
    
    # Combine datasets
    #density_df = pd.concat([df_greenland, df_antarctica], ignore_index=True) #maybe useful for 'both' option later
    density_df = df.copy()
    
    print(f"Loaded {len(df):,} {region} density measurements")
    
    # Load profile names
    profile_names_file = os.path.join(data_path, 'SUMup_2024_density_profile_names.tsv')

    if os.path.exists(profile_names_file):
        profile_names_df = pd.read_csv(profile_names_file, sep='\t').set_index('key')
        print(f"Loaded {len(profile_names_df)} profile names")
    else:
        print("Warning: Profile names file not found")
        profile_names_df = pd.DataFrame()
    
    # Load references
    references_file = os.path.join(data_path, 'SUMup_2024_density_references.tsv')
    if os.path.exists(references_file):
        references_df = pd.read_csv(references_file, sep='\t').set_index('key')
        print(f"Loaded {len(references_df)} references")
    else:
        print("Warning: References file not found")
        references_df = pd.DataFrame()
    
    return density_df, profile_names_df, references_df

def def_models():
        
    def logarithmic(x, a, b, c):
        return a * np.log(x + 1) + b * x + c

    def log_minimization(x, a, b, c, target):
        return a * np.log(x + 1) + b * x + c - target

    def power_law(x, a, b, c):
        return a * x**b + c

    def power_law_minimization(x, a, b, c, target):
        return a * x**b + c - target

    def linear(x, m, b):
        return m * x + b
    
    def linear_minimization(x, m, b, target):
        return m * x + b - target


    models = [
        (power_law, "Power law", "a*x^b + c", 'blue', power_law_minimization),
        (logarithmic, "Logarithmic", "a*log(x+1) + b*x + c", 'green', log_minimization),
        (linear, "Linear", "m*x + b", 'orange', linear_minimization)
    ]
    
    return models

def calculate_depth_to_density_sumup(depth, density, core_name, depth_dependence=True, do_plot=False):
    """
    Interpolate to find depth at target density
    
    Parameters:
    depth (array): Depth measurements (should be midpoint values)
    density (array): Corresponding density measurements
    core_name (str): name of core
    depth_dependence (bool): Whether to include depth dependence (only use {550: min density- 600 kg/m3; 830: 550 kg/m3 - max density} when fitting
    do_plot (bool): Whether to plot the fit
    
    Returns:
    depths_at_target_densities (array of floats or nans): Depths at 550 and 830 kg/m³
    failed (array of ints): 0 if successful, 1 if failed
    r2 (array of floats): R² value of the fit
    
    """

    print(f"Processing {core_name}....\n")
    
    target_densities = [550, 830]
    colors = ['blue', 'red']
    depths_at_target_densities = [np.nan,np.nan]
    failed = [0,0]
    r2 = [np.nan, np.nan]
    saved = 0

    models=def_models()

    # Remove NaN values
    valid_mask = ~(np.isnan(depth) | np.isnan(density))

    if valid_mask.sum() < 10:
        
        failed = [1,1]
        print(" - Too few valid points ("+str(valid_mask.sum())+")")
        return depths_at_target_densities, failed, r2
            
    depth = np.array(depth[valid_mask])
    density = np.array(density[valid_mask])
    
    # Sort by depth for interpolation
    sort_idx = np.argsort(depth)
    depth_sorted = depth[sort_idx]
    density_sorted = density[sort_idx]
    
    # Check if target density is within the range of measured densities
    min_density = density_sorted.min()
    max_density = density_sorted.max()

    #Establishes figure
    if do_plot:
        
        fig,ax = plt.subplots(figsize=(5,5))
        ax.plot(density, depth, 'o', label='Core', color='grey')
        plt.gca().invert_yaxis()
        ax.set_title(f'{core_name}')

    # target density out of range
    for idx, rho_i in enumerate(target_densities):

        # check if target density is within the range of measured densities and/or has more than 10 observatinos
        if rho_i < min_density or rho_i > max_density:
            failed[idx] = 1
            print("  - Target density out of range (" + str(rho_i) + " not in [" + str(min_density) + ", " + str(max_density) + "])")
        
        if all(failed):
            return depths_at_target_densities, failed, r2

    # Try fit to a logarithmic or power law curve.
    func, name, expression, color, func_min = models[1] # only log
    
    for idx, rho_i in enumerate(target_densities):

        # if this density is not in the data, skip
        if failed[idx] == 1:
            continue

        try: 

            # only uses specific depths to fit log function -- improves quality of depth estimates
            if depth_dependence:
                if rho_i == 550:
                    fit_mask = density_sorted <= 600  # fit only to densities up to 600
                else:
                    fit_mask = density_sorted >= 550  # fit only to densities above 600
            else:
                fit_mask = np.ones_like(density_sorted, dtype=bool)  # use all data

            if len(depth_sorted[fit_mask]) < 10:
                raise Exception("Too few points to fit for {rho_i}")

            # Fit the model
            popt, pcov = curve_fit(func, depth_sorted[fit_mask], density_sorted[fit_mask], maxfev=10000)

            # Calculate fitted values
            depth_smooth = np.linspace(depth_sorted.min(), depth_sorted.max(), 200)
            density_fit = func(depth_sorted, *popt)
            density_smooth = func(depth_smooth, *popt)

            # Calculate R²
            residuals = density_sorted[fit_mask] - density_fit[fit_mask]
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((density_sorted[fit_mask] - np.mean(density_sorted[fit_mask]))**2)
            r2[idx] = round(1 - (ss_res / ss_tot), 2)

            failed[idx] = 0

            if r2[idx]<0.85:
                failed[idx] = 1
                raise Exception("Low r2; failed")
        
        except Exception as e:
            failed[idx] = 1
            print(f"Model fitting failed for {rho_i}")
            continue
        

        # try to find the depth at target density using the fit model
        try:
            depth_at_rho = brentq(
                func_min,
                depth_sorted[fit_mask].min(),
                depth_sorted[fit_mask].max(),
                args=(*popt, rho_i)
            )

            depths_at_target_densities[idx] = round(depth_at_rho, 2)
            
        except Exception as e:
            print(f"  - Root finding failed for {rho_i}.")
            failed[idx] = 1
            pass
            
        if do_plot & (not failed[idx]):
            
            ax.plot(density_smooth, depth_smooth, '-', label=f'{rho_i} model, r2={r2[idx]}', color=colors[idx])
            ax.plot(rho_i, depth_at_rho, 'x', color=colors[idx], markersize=10)
            ax.plot(np.ones(len(depth_smooth))*rho_i, depth_smooth, '--', color='lightgray')

            ax.legend()
            ax.set_xlabel('Density (kg/m³)')
            ax.set_ylabel('Depth (m)')
            ax.set_xlim([250,950])

            if depth_dependence:
                saved = 1
                plt.savefig(f"../figures/sumup-core-fits/dd/{core_name}-dd.png")
            else:
                saved = 1
                plt.savefig(f"../figures/sumup-core-fits/whole/{core_name}-whole.png")
                
    if (saved==0) & do_plot:
        plt.savefig(f"../figures/sumup-core-fits/just-observations/{core_name}.png")
        
    plt.close('all')
    return depths_at_target_densities, failed, r2

def process_cores_sumup(density_df, profile_names_df, references_df, depth_dependence=True, do_plot=False):
    """
    Process all cores to find depths at target densities
    
    Parameters:
    density_df (DataFrame): Combined density measurements
    profile_names_df (DataFrame): Profile names lookup
    references_df (DataFrame): References lookup
    
    Returns:
    DataFrame: Core information with depths at 550 and 830 kg/m³
    """
    
    print("Processing individual cores...")
    
    results = []
    
    # Group by profile_key to process each core separately
    grouped = density_df.groupby('profile_key')
    
    total_profiles = len(grouped)
    processed = 0
    
    for profile_key, group in grouped:
        processed += 1
        if processed % 100 == 0:
            print(f"Processed {processed}/{total_profiles} profiles...")
        
        # Get core metadata
        profile_name = profile_names_df.loc[profile_key, 'profile'] if profile_key in profile_names_df.index else f"Profile_{profile_key}"
        
        # Get reference info
        reference_key = group['reference_key'].iloc[0]
        reference_short = references_df.loc[reference_key, 'reference_short'] if reference_key in references_df.index else f"Ref_{reference_key}"
        
        # Get basic info
        date = group['timestamp'].iloc[0] if pd.notna(group['timestamp'].iloc[0]) else "Unknown"
        region = group['region'].iloc[0]
        latitude = group['latitude'].iloc[0]
        longitude = group['longitude'].iloc[0]
        elevation = group['elevation'].iloc[0]
        
        # Use midpoint for depth if available, otherwise calculate from start/stop
        if 'midpoint' in group.columns and group['midpoint'].notna().any():
            depths = group['midpoint'].values
        else:
            # Calculate midpoint from start_depth and stop_depth if available
            if 'start_depth' in group.columns and 'stop_depth' in group.columns:
                depths = (group['start_depth'] + group['stop_depth']) / 2
            else:
                print(f"Warning: No depth information available for profile {profile_key}")
                continue
        
        densities = group['density'].values
        
        if profile_name is np.nan:
            core_name = "Key: " + str(profile_key)
        else:
            core_name = profile_name
            
        # Fit a model to find depths at target densities
        depth_at_density, failed, r2 = calculate_depth_to_density_sumup(depths, densities, core_name, depth_dependence=depth_dependence, do_plot=do_plot)
        
        result_row = {
            'core_name': str(profile_name),
            'profile_key': profile_key,
            'date': date,
            'region': region,
            'latitude': latitude,
            'longitude': longitude,
            'elevation': elevation,
            'citation': reference_short,
            'depth_to_550': depth_at_density[0],
            'depth_to_830': depth_at_density[1],
            'failed_fit_550': failed[0],
            'failed_fit_830': failed[1],
            'r2_550': r2[0],
            'r2_830': r2[1]
        }
        
        results.append(result_row)
    
    print(f"Processed {processed} profiles total")
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Add year column to better match with the metadata from Peter's standardized cores
    results_df["date"] = pd.to_datetime(results_df["date"])
    results_df["year"] = results_df["date"].dt.year

    # Add summary statistics
    print("\nSummary:")
    print(f"Total cores processed: {len(results_df)}")
    print(f"Cores reaching 550 kg/m³: {results_df['depth_to_550'].notna().sum()}")
    print(f"Cores reaching 830 kg/m³: {results_df['depth_to_830'].notna().sum()}")
    print(f"Cores reaching both densities: {(results_df['depth_to_550'].notna() & results_df['depth_to_830'].notna()).sum()}")
    
    return results_df

def save_results(results_df, output_path, output_prefix='sumup_density_depths'):
    """
    Save results to separate CSV files for each region
    Only includes cores that reach at least 550 kg/m³ density
    
    Parameters:
    results_df (DataFrame): Results to save
    output_prefix (str): Prefix for output filenames
    """
    
    # Filter to only include cores with either depth
    filtered_df = results_df[results_df['depth_to_830'].notna()+results_df['depth_to_550'].notna()].copy()
    
    print(f"\nFiltering results:")
    print(f"Total cores processed: {len(results_df)}")
    print(f"Cores with at least one depth: {len(filtered_df)}")
    
    # Split filtered data by region
    
    #maybe useful for 'both' option later
    #greenland_df = filtered_df[filtered_df['region'] == 'Greenland'].copy()
    #antarctica_df = filtered_df[filtered_df['region'] == 'Antarctica'].copy()
    
    # Sort by core name within each region
    #maybe useful for 'both' option later``
    #greenland_df = greenland_df.sort_values('core_name')
    #antarctica_df = antarctica_df.sort_values('core_name')
    filtered_df = filtered_df.sort_values('core_name')
    
    # Save separate files
    try:
        region = filtered_df['region'].unique()[0]
    except:
        raise ValueError(f"Either no region specified or multiple regions found in filtered data. Regions: {region}")

    file_path = f"{output_path}{output_prefix}_{region}.csv"
    filtered_df.to_csv(file_path, index=False, float_format='%.3f')

    print(f"{region} cores: {len(filtered_df)}")
    print(f"  - Reaching 550 kg/m³: {filtered_df['depth_to_550'].notna().sum()}")
    print(f"  - Reaching 830 kg/m³: {filtered_df['depth_to_830'].notna().sum()}")
    # Save regional files
    
    # Print some example results for each region
    display_cols = ['core_name', 'date', 'depth_to_550', 'depth_to_830']
    
    print(f"\nSample of processed cores:")
    print(filtered_df[display_cols].head().to_string(index=False))

def process_sumup_main(data_path = "../data/sumup/source/", output_path = "../data/sumup/processed/", domain='greenland', do_plot=False):
    """
    Main function
    
    Inputs:
        - data_path (str): Path to folder containing SUMup CSV files
        - output_path (str): Path to save output files
        - domain (str): 'greenland' or 'antarctica'

    Outputs:
        - (none) saves processed CSV files to output_path
    """
    
    # ==========================================
    # CONFIGURATION - UPDATE THIS PATH
    # ==========================================

    output_prefix = "sumup_550_830_density_depths" # Will create _greenland.csv, _antarctica.csv, and _combined.csv files
    
    # Check if data path exists
    if not os.path.exists(data_path):
        print(f"Data path does not exist: {data_path}")
        download_sumup_data_instructions()
        return
    
    try:
        # Load data
        density_df, profile_names_df, references_df = load_sumup_data(data_path, domain)
        
        # Process cores
        results_df = process_cores_sumup(density_df, profile_names_df, references_df, depth_dependence=True, do_plot=do_plot)
        results_df["source"] = "SUMUP 2024"
        results_df["source_doi"] = "doi:10.18739/A2M61BR5M"

        # now we want to do some cleanup -- attempt to refit the profiles that have no 550 or 830 values
        fit_again_names = ["ACT10A", "ACT11C", "ACT11D", "Dye 2 (Raven) Core A", "Site J"]
        fit_again_keys_sel = [562] + results_df[results_df["core_name"].isin(fit_again_names)].profile_key.to_list()
        fit_again_keys = fit_again_keys_sel + results_df[(results_df["depth_to_550"].isna()) | (results_df["depth_to_830"].isna())].profile_key.to_list()

        fit_again_density = density_df[density_df["profile_key"].isin(fit_again_keys)]
        fit_again_results_df = process_cores_sumup(fit_again_density, profile_names_df, references_df, depth_dependence=False, do_plot=do_plot)

        # only fills in new fit_again (whole profile) depths where whole density fit r2 > 0.8 AND no depth-dependent value was found
        new_550_df = fit_again_results_df[~fit_again_results_df["depth_to_550"].isna()].copy()
        new_830_df = fit_again_results_df[~fit_again_results_df["depth_to_830"].isna()].copy()

        for idx, row in results_df.iterrows():

            if np.isin(row["profile_key"],new_550_df["profile_key"]) & np.isnan(row["depth_to_550"]):

                new_depth = new_550_df.depth_to_550[new_550_df["profile_key"]==row["profile_key"]].values[0]
                results_df.loc[idx,"depth_to_550"] = row.depth_to_550
            
            if (np.isin(row["profile_key"],new_830_df["profile_key"])) & np.isnan(row["depth_to_830"]):
                
                new_depth = new_830_df.depth_to_830[new_830_df["profile_key"]==row["profile_key"]].values[0]
                results_df.loc[idx,"depth_to_830"] = row.depth_to_830

        # # also pull in data from the new fits for hand selected cores that were fit with the depth-dependence scheme
        good_whole_fits = ["ACT10A", "ACT10B", "ACT11A", "ACT11C", "ACT11D", "Dye 2 (Raven) Core A"]
        good_whole_keys = [562] + results_df[results_df["core_name"].isin(fit_again_names)].profile_key.to_list()

        for idx, row in results_df.iterrows():

            # if one of the selected cores
            if np.isin(row['profile_key'],good_whole_keys):

                if np.isin(row['profile_key'], new_550_df['profile_key']):
                    new_depth_550 = new_550_df.depth_to_550[new_550_df['profile_key']==row['profile_key']].values[0]
                    if ~np.isnan(new_depth_550):
                        results_df.loc[idx,"depth_to_550"] = new_depth_550

                if np.isin(row['profile_key'], new_830_df['profile_key']):
                    new_depth_830 = new_830_df.depth_to_830[new_830_df['profile_key']==row['profile_key']].values[0]
                    if ~np.isnan(new_depth_830):
                        results_df.loc[idx,"depth_to_830"] = new_depth_830

        # drop hand-selected cores because bad data slips through

        drop_names = ["FS2_12m"]
        drop_keys = [1951,955] + results_df[results_df["core_name"].isin(drop_names)].profile_key.to_list()
        drop_idxs = results_df[results_df["profile_key"].isin(drop_keys)].index.values
        results_df = results_df.drop(index = drop_idxs).reset_index(drop=True)

        # Save results
        save_results(results_df, output_path, output_prefix)
        
        print(f"\n✓ Processing complete!")

    except Exception as e:

        print(f"Error: {e}")
        print("\nMake sure you have downloaded the SUMup data files.")
        download_sumup_data_instructions()

if __name__ == "__main__":
    process_sumup_main()