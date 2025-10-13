import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import os
import warnings
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

def interpolate_depth_to_density(depth_values, density_values, target_density):
    """
    Interpolate to find depth at target density
    
    Parameters:
    depth_values (array): Depth measurements (should be midpoint values)
    density_values (array): Corresponding density measurements
    target_density (float): Target density to find depth for (kg/m³)
    
    Returns:
    float or NaN: Depth at target density, or NaN if interpolation not possible
    """
    
    # Remove NaN values
    valid_mask = ~(np.isnan(depth_values) | np.isnan(density_values))
    if valid_mask.sum() < 2:
        return np.nan
        
    depth_clean = depth_values[valid_mask]
    density_clean = density_values[valid_mask]
    
    # Sort by depth for interpolation
    sort_idx = np.argsort(depth_clean)
    depth_sorted = depth_clean[sort_idx]
    density_sorted = density_clean[sort_idx]
    
    # Check if target density is within the range of measured densities
    min_density = density_sorted.min()
    max_density = density_sorted.max()
    
    if target_density < min_density or target_density > max_density:
        return np.nan
    
    # Check if densities are monotonically increasing with depth
    # (required for meaningful interpolation)
    if not np.all(np.diff(density_sorted) >= 0):
        # Try to handle non-monotonic data by using only the deepest measurement
        # for each density (common in firn cores)
        unique_densities = []
        unique_depths = []
        for i, dens in enumerate(density_sorted):
            if dens not in unique_densities:
                unique_densities.append(dens)
                unique_depths.append(depth_sorted[i])
        
        if len(unique_densities) < 2:
            return np.nan
            
        density_sorted = np.array(unique_densities)
        depth_sorted = np.array(unique_depths)
    
    try:
        # Linear interpolation
        f = interp1d(density_sorted, depth_sorted, kind='linear', 
                    bounds_error=False, fill_value=np.nan)
        return float(f(target_density))
    except:
        return np.nan

def process_cores(density_df, profile_names_df, references_df):
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
            'citation': reference_short,
            'n_measurements': len(group),
            'max_depth': np.nanmax(depths),
            'max_density': np.nanmax(densities),
        }
        
        # Interpolate to find depths at target densities
        for target_density in target_densities:
            depth_at_density = interpolate_depth_to_density(depths, densities, target_density)
            result_row[f'depth_to_{target_density}'] = depth_at_density
        
        results.append(result_row)
    
    print(f"Processed {processed} profiles total")
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Add year column to better match with the metadata from Peter's standardized cores
    results_df["date"] = pd.to_datetime(results_df["date"])
    results_df["year"] = results_df["date"].dt.year
    results_df["source"] = "SUMUP 2024"
    results_df["source_doi"] = "doi:10.18739/A2M61BR5M"

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

def process_sumup_main():
    """Main function"""
    
    # ==========================================
    # CONFIGURATION - UPDATE THIS PATH
    # ==========================================
    data_path = "../data/sumup/raw/"  # Update this to your SUMup data folder
    output_path = "../data/sumup/processed/"  # Path to save output files (can be same as data_path)
    output_prefix = "sumup_550_830_density_depths"  # Will create _greenland.csv, _antarctica.csv, and _combined.csv files
    domain = "greenland"
    
    # Check if data path exists
    if not os.path.exists(data_path):
        print(f"Data path does not exist: {data_path}")
        download_sumup_data_instructions()
        return
    
    try:
        # Load data
        density_df, profile_names_df, references_df = load_sumup_data(data_path, domain)
        
        # Process cores
        results_df = process_cores(density_df, profile_names_df, references_df)
        
        # Save results
        save_results(results_df, output_path, output_prefix)
        
        print(f"\n✓ Processing complete!")

    except Exception as e:

        print(f"Error: {e}")
        print("\nMake sure you have downloaded the SUMup data files.")
        download_data_instructions()

if __name__ == "__main__":
    process_sumup_main()