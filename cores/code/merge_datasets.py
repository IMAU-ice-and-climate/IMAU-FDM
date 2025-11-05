import pandas as pd
import numpy as np
import os
from datetime import datetime

def load_sumup_dataset(sumup_file):
    """Load pre-processed SUMUP dataset"""
    print(f"\nLoading SUMUP dataset from {sumup_file}")
    
    sumup_df = pd.read_csv(sumup_file)

    print(f"Loaded {len(sumup_df)} SUMUP records")
    print("SUMUP columns:", sumup_df.columns.tolist())

    
    return sumup_df

def load_standardized_dataset(standardized_file):
    """Load pre-processed standardized dataset"""
    print(f"\nLoading standardized dataset from {standardized_file}")
    
    std_df = pd.read_csv(standardized_file)

    print(f"Loaded {len(std_df)} standardized records")
    print("Standardized columns:", std_df.columns.tolist())

    
    return std_df

def identify_duplicates(std_df, sumup_df, location_tolerance=0.1, year_tolerance=0):
    """Identify potential duplicates between new cores and SUMUP data"""
    print(f"\nIdentifying duplicates...")
    
    duplicates = []

    sumup_cores_no_dupes = sumup_df.copy(deep=True) # default to drop dupes from sumup dataset and keep sumup dataset, since std was handpicked
        
    for row_i, new_core in sumup_df.iterrows():
        
        # Check for name matches
        name_matches = std_df[std_df['core_name'].str.lower() == str(new_core['core_name']).lower()]
        
        if name_matches.empty:
            # Check for location matches
            if not pd.isna(new_core['latitude']) and not pd.isna(new_core['longitude']):
                lat_diff = abs(std_df['latitude'] - new_core['latitude'])
                lon_diff = abs(std_df['longitude'] - new_core['longitude'])
                location_matches = std_df[(lat_diff <= location_tolerance) & (lon_diff <= location_tolerance)]
                
                # Further filter by year if available
                if not pd.isna(new_core['year']):
                    year_diff = abs(location_matches['year'] - new_core['year'])
                    location_matches = location_matches[
                        pd.isna(year_diff) | (year_diff <= year_tolerance)
                    ]
            else:
                location_matches = pd.DataFrame()
        
        all_matches = pd.concat([name_matches, location_matches])
        
        if not all_matches.empty:
            for _, match in all_matches.iterrows():
                duplicate_info = {
                    'new_core_name': new_core['core_name'],
                    'new_core_year': new_core['year'],
                    'new_core_lat': new_core['latitude'],
                    'new_core_lon': new_core['longitude'],
                    'new_core_depth_to_550': new_core.get('depth_to_550', np.nan),
                    'new_core_ref': new_core.get('citation', ''),
                    'sumup_name': match['core_name'],
                    
                    'sumup_year': match.get('year', np.nan),
                    'sumup_lat': match['latitude'],
                    'sumup_lon': match['longitude'],
                    'sumup_depth_to_550': match.get('depth_to_550', np.nan),
                    'sumup_ref': match.get('citation', ''),
                    'match_type': 'core_name' if str(new_core['core_name']).lower() == str(match['core_name']).lower() else 'location'

                }
                duplicates.append(duplicate_info)

            sumup_cores_no_dupes = sumup_cores_no_dupes.drop(index=row_i) # drops duplicates from the sumup core set

    duplicates_df = pd.DataFrame(duplicates)
    print(f"Found {len(duplicates_df)} potential duplicate pairs")

    
    return duplicates_df, sumup_cores_no_dupes

def create_merged_dataset(std_df, sumup_df):
    """Create final combined dataset in SUMUP format"""
    print(f"\nCreating final dataset...")
    
    # Ensure both datasets have the same columns as SUMUP format
    column_names = ['core_name', 'latitude', 'longitude', 'elevation', 'year', 'depth_to_550', 'depth_to_830', 'r2_550', 'r2_830', 'failed_550','failed_830','source', 'citation', 'region']
    
    # Select and reorder columns
    std_final = std_df[column_names].copy(deep=True)
    std_final['r2'] = np.nan  # Placeholder for r2 in standardized data
    
    # Prepare SUMUP data
    sumup_final = sumup_df[column_names].copy(deep=True)
    
    # Combine datasets
    final_df = pd.concat([std_final, sumup_final], ignore_index=True)
    
    # Remove cores missing depth_to_550
    initial_count = len(final_df)
    final_df = final_df.dropna(subset=['depth_to_550'])
    removed_count = initial_count - len(final_df)
    
    print(f"Removed {removed_count} cores missing depth_to_550")
    print(f"Final dataset contains {len(final_df)} cores")
    
    # Sort by study then name
    final_df = final_df.sort_values(['year','core_name']).reset_index(drop=True)
    
    return final_df

def merge_datasets(do_run_sumup_processing=False, do_run_std_processing=False, drop_duplicates=True, return_df=True, data_dir="../data/"):

    """
        Merges the sumup and standardized datasets; more can be added as long as core file structure stays the same.

        # INPUTS set pre-processing options
        do_run_sumup_processing = False  # Set to True to re-process SUMUP data from raw files
        do_run_std_processing = False  # Set to True to process standardized core data

        drop_duplicates = True  # Whether to drop duplicates from standardized dataset found in SUMUP dataset when merging 
    """

    # set filepaths
    
    region = "greenland"
    current_year = str(datetime.now().year)

    sumup_file = f"{data_dir}sumup/processed/sumup_550_830_density_depths_{region}.csv"
    std_file = f"{data_dir}PKM/processed/PKM_550_830_density_depths_{region}.csv"
    
    merged_file = f"{data_dir}merged/MERGED_CORE_LIST_{region}_{current_year}.csv"

    #re-run processing if needed

    if do_run_sumup_processing:
        print("\nRe-processing SUMUP data from raw files...")

        try:
            import process_sumup_cores
            process_sumup_cores.process_sumup_main()  # This will re-generate the SUMUP dataset file
        except Exception as e:
            print(f"Error re-processing SUMUP data: {e}")
            return
        
    if do_run_std_processing:
        print("\nRe-processing standardized data from raw files...")

        try:
            import process_standardized_cores
            process_standardized_cores.process_std_main()  # This will re-generate the standardized dataset file
        except Exception as e:
            print(f"Error re-processing standardized data: {e}")
            return
        
    # Load existing datasets
            
    if os.path.exists(sumup_file):
        try:
            print(f"Loading sumup data from file... {sumup_file}")
            sumup_df = load_sumup_dataset(sumup_file)
        except Exception as e:
            print(f"Error loading SUMUP data: {e}")
            return
        
    if os.path.exists(std_file):
        try:
            print(f"Loading sumup data from file... {std_file}")
            std_df = load_standardized_dataset(std_file)
        except Exception as e:
            print(f"Error loading standardized data: {e}")
            return
        
    # Identify duplicates
    duplicates_df, std_cores_no_dupes = identify_duplicates(std_df, sumup_df)

    if drop_duplicates:
        final_df = create_merged_dataset(std_cores_no_dupes, sumup_df)
    else:
        final_df = create_merged_dataset(std_df, sumup_df)

    # Save out final merged dataset
    try:
        final_df.to_csv(merged_file, index=False)
        print(f"\n=== MERGING COMPLETE ===")
        print(f"Final merged dataset saved to: {merged_file}")
        print(f"Cores with 550 kg/m³ horizon: {len(final_df[~final_df['depth_to_550'].isna()])}")
        print(f"Cores with 830 kg/m³ horizon: {len(final_df[~final_df['depth_to_830'].isna()])}")
    except Exception as e:
        print(f"Error saving final merged dataset: {e}")

    if return_df:
        return final_df, duplicates_df, sumup_df, std_df

if __name__ == "__main__":
    merge_datasets()
