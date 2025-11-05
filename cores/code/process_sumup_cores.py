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
https://arcticdata.io/catalog/view/doi:10.18739/A2M61BR5M

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


    models = [
        (power_law, "Power law", "a*x^b + c", 'blue', power_law_minimization),
        (logarithmic, "Logarithmic", "a*log(x+1) + b*x + c", 'green', log_minimization)
    ]
    
    return models

def interpolate_depth_to_density_sumup(depth, density, target_density, name, do_plot=False):
    """
    Interpolate to find depth at target density
    
    Parameters:
    depth (array): Depth measurements (should be midpoint values)
    density (array): Corresponding density measurements
    target_density (float): Target density to find depth for (kg/m³)
    
    Returns:
    depth_at_rho: float or NaN: Depth at target density, or NaN if interpolation not possible
    

    """

    models=def_models()

    failed = 0
    depth_at_rho=np.nan

    best_r2 = -np.inf
    best_model = None


    # Remove NaN values
    valid_mask = ~(np.isnan(depth) | np.isnan(density))

    if valid_mask.sum() < 2:
        depth_at_rho = np.nan
        failed = 1
        return depth_at_rho, failed, best_r2
            
    depth_clean = depth[valid_mask]
    density_clean = density[valid_mask]
    
    # Sort by depth for interpolation
    sort_idx = np.argsort(depth_clean)
    depth_sorted = depth_clean[sort_idx]
    density_sorted = density_clean[sort_idx]
    
    # Check if target density is within the range of measured densities
    min_density = density_sorted.min()
    max_density = density_sorted.max()

    if target_density < min_density or target_density > max_density or len(depth_sorted)<=10:
        
        #print("  - Target density out of range")
        failed = 1
        return depth_at_rho, failed, best_r2

        
    # Try fit to a logarithmic or power law curve.

    for idx, (func, name, expression, color, func_min) in enumerate(models):

        try: 

            # Fit the model
            popt, pcov = curve_fit(func, depth_sorted, density_sorted)

            # Calculate fitted values
            depth_smooth = np.linspace(depth_sorted.min(), depth_sorted.max(), 200)
            density_fit = func(depth_sorted, *popt)
            density_smooth = func(depth_smooth, *popt)

            # Calculate R²
            residuals = density_sorted - density_fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((density_sorted - np.mean(density_sorted))**2)
            r2 = 1 - (ss_res / ss_tot)

            failed = 0
            
            if r2 > best_r2:
                
                best_r2 = r2
                best_model = (func, popt, name, expression, round(best_r2,3), func_min)                
        
        except Exception as e:
            pass
            #print(f"Model fitting failed for {name}: {e}")
            

    if best_r2 < 0.85:
        failed = 1
        #print(f"  R² ({best_r2:.3f}) too low.")
        return depth_at_rho, failed, round(best_r2, 2)
        
    
    # try to find the depth at target density using the best model
    try:
        depth_at_rho = brentq(
            best_model[5],
            depth_sorted.min(),
            depth_sorted.max(),
            args=(*best_model[1], target_density)
        )

        
        #print("fit succeeded")
        
    except Exception as e:
        #print(f"  - Root finding failed for best model: {best_model[2]}.")
        failed = 1
            
    if do_plot:

        fig,ax = plt.subplots(figsize=(5,5))
        ax.plot(density, depth, 'o', label='Core', color='grey')
        
        if failed == 0:
            ax.plot(density_smooth, depth_smooth, '-', label='Log fit', color='blue')
        
        ax.plot(target_density, depth_at_rho, 'rx', label=f'Depth at {target_density} kg/m³', markersize=10)

        plt.legend()
        plt.gca().invert_yaxis()
        plt.xlabel('Density (kg/m³)')
        plt.ylabel('Depth (m)')
        plt.title(f'Depth at {target_density} kg/m³: {depth_at_rho:.2f} m')
        plt.show()
        plt.savefig(f"/home/nld4814/perm/cores/sumup/data/processed/figs/{name}_{str(target_density)}.png")

    
    return round(depth_at_rho, 2), failed, round(best_r2, 2)

def process_cores_sumup(density_df, profile_names_df, references_df, do_plot=False):
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
    
    target_densities = [550, 830]  # kg/m³
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
        
        # Find depths for target densities
        result_row = {
            'core_name': profile_name,
            'profile_key': profile_key,
            'date': date,
            'region': region,
            'latitude': latitude,
            'longitude': longitude,
            'elevation': elevation,
            'citation': reference_short
        }
        
        # Interpolate to find depths at target densities
        for target_density in target_densities:

            depth_at_density, failed, r2 = interpolate_depth_to_density_sumup(depths, densities, target_density, profile_name, do_plot=do_plot)

            result_row[f'depth_to_{target_density}'] = depth_at_density
            result_row[f'failed_fit_{target_density}'] = failed
            result_row[f'r2_{target_density}'] = r2
        
        results.append(result_row)
    
    print(f"Processed {processed} profiles total")
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Add year column to better match with the metadata from Peter's standardized cores
    results_df["date"] = pd.to_datetime(results_df["date"])
    results_df["year"] = results_df["date"].dt.year
    results_df["source"] = "SUMUP 2024"
    results_df["source_doi"] = "doi:10.18739/A2M61BR5M"

    results_df = results_df[results_df['failed_fit_550'] == 0]

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
    
    # Filter to only include cores that reach at least 550 kg/m³
    filtered_df = results_df[results_df['depth_to_550'].notna()].copy()
    
    print(f"\nFiltering results:")
    print(f"Total cores processed: {len(results_df)}")
    print(f"Cores reaching 550 kg/m³: {len(filtered_df)}")
    print(f"Cores filtered out: {len(results_df) - len(filtered_df)}")
    
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

def process_sumup_main(data_path = "../data/sumup/source/", output_path = "../data/sumup/processed/", domain='greenland'):
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
    #data_path = "../data/sumup/raw/"  # Update this to your SUMup data folder
    #output_path = "../data/sumup/processed/"  # Path to save output files (can be same as data_path)
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
        results_df = process_cores_sumup(density_df, profile_names_df, references_df)
        
        # Save results
        save_results(results_df, output_path, output_prefix)
        
        print(f"\n✓ Processing complete!")

    except Exception as e:

        print(f"Error: {e}")
        print("\nMake sure you have downloaded the SUMup data files.")
        download_sumup_data_instructions()

if __name__ == "__main__":
    process_sumup_main()