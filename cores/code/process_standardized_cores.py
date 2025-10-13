import pandas as pd
import numpy as np
import os
import glob
import warnings
warnings.filterwarnings('ignore')

"""
Process std.txt files from standardized core data (originally processed by Peter Kuipers-Munneke) to create a combined dataset in a standardized format
"""

# STD.TXT
def load_std_metadata(metadata_file):

    """Load std.txt firn core metadata from Excel file"""
    print(f"Loading metadata from {metadata_file}")
    
    # Load the '62 cores' sheet (index 1) which contains Greenland data
    df = pd.read_csv(metadata_file)
    
    # Clean up column names and data
    df = df.dropna(subset=['Name'])  # Remove rows without names
    df = df.reset_index(drop=True)
    
    print(f"Loaded metadata for {len(df)} cores")
    print("Sample core names:", df['Name'].head().tolist())
    
    return df

def calculate_depth_to_density(depths, densities, target_density):
    """Calculate depth where density first reaches target value using interpolation"""
    #TODO: update this interpolation function to match the sumup one
    
    depths = np.array(depths)
    densities = np.array(densities)
    
    # Sort by depth
    sort_idx = np.argsort(depths)
    depths_sorted = depths[sort_idx]
    densities_sorted = densities[sort_idx]
    
    # Check if target density is reached
    if np.max(densities_sorted) < target_density:
        return np.nan
    
    # Find first point where target density is reached or exceeded
    idx = np.where(densities_sorted >= target_density)[0]
    if len(idx) == 0:
        return np.nan
    
    first_idx = idx[0]
    
    # If first measurement already exceeds target
    if first_idx == 0:
        return depths_sorted[0]
    
    # Linear interpolation between the two points
    prev_depth = depths_sorted[first_idx - 1]
    prev_density = densities_sorted[first_idx - 1]
    curr_depth = depths_sorted[first_idx]
    curr_density = densities_sorted[first_idx]
    
    # Interpolate
    if curr_density == prev_density:
        return prev_depth
    
    interpolated_depth = prev_depth + (target_density - prev_density) * (curr_depth - prev_depth) / (curr_density - prev_density)
    
    return round(interpolated_depth, 2)

def import_and_process_std_text_files(text_files_pattern, metadata_df):
    """Process all firn core std text files and calculate depth to 550 and 830"""
    print(f"\nProcessing text files matching pattern: {text_files_pattern}")
    
    text_files = glob.glob(text_files_pattern)
    print(f"Found {len(text_files)} text files")
    
    processed_cores = []
    
    for file_path in text_files:
        # Extract core name from filename
        filename = os.path.basename(file_path)
        core_name = filename.replace('.std.txt', '').replace('.txt', '')
        
        try:
            # Read the text file
            data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['depth', 'density'])
            
            # Find matching metadata
            metadata_match = metadata_df[metadata_df['Name'].astype(str).str.lower() == core_name.lower()]
            
            if metadata_match.empty:
                # Try numeric match for cases like '6345'
                try:
                    core_num = float(core_name)
                    metadata_match = metadata_df[metadata_df['Name'] == core_num]
                except:
                    pass
            
            if metadata_match.empty:
                print(f"Warning: No metadata found for core {core_name}")
                continue
            
            metadata_row = metadata_match.iloc[0]
            
            # Calculate depth to density horizons
            depth_550 = calculate_depth_to_density(data['depth'], data['density'], 550)
            depth_830 = calculate_depth_to_density(data['depth'], data['density'], 830)
            
            # Create core record
            core_record = {
                'core_name': str(metadata_row['Name']),
                'latitude': metadata_row['Latitude'],
                'longitude': metadata_row['Longitude'],
                'elevation': metadata_row.get('Elevation', np.nan),
                'year': metadata_row.get('Year', np.nan),
                'depth_to_550': depth_550,
                'depth_to_830': depth_830,
                'source': '"PKM standardized cores"',
                'method': metadata_row.get('Method', ''),
                'citation': metadata_row.get('Citation', ''),
                'measurement_count': len(data),
                'source_file': filename,
                "region": "Greenland"
            }
            
            processed_cores.append(core_record)
            print(f"Processed {core_name}: {len(data)} measurements, depth_550={depth_550}, depth_830={depth_830}")
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    cores_df = pd.DataFrame(processed_cores)
    
    # fix elevation column, which includes some elevations of "X"
    cores_df['elevation'] = cores_df['elevation'].replace('X', np.nan)
    cores_df["elevation"] = cores_df["elevation"].replace('NaN', np.nan)
    cores_df["elevation"] = cores_df["elevation"].astype(float)

    print(f"\nProcessed {len(cores_df)} cores successfully")
    
    return cores_df

def process_std_main():
    """Main processing function"""
    
    domain = "greenland"
    # File paths - UPDATE THESE PATHS AS NEEDED
    std_core_dir = "/../data/cores/raw/PKM/"

    metadata_file = std_core_dir+"cleaned-std-metadata.csv"
    text_files_pattern = std_core_dir+"standardized/*.txt"  # or "*.std.txt" or "/path/to/cores/*.txt"
    output_file_std = f"{std_core_dir}processed/PKM_550_830_density_depths_{domain}.csv"
    
    print("=== Processing std.txt files from PKM's dataset ===")
    
    # Step 1: Load metadata
    try:
        metadata_df = load_std_metadata(metadata_file)
    except Exception as e:
        print(f"Error loading metadata: {e}")
        return
    
    # Step 2: Process text files
    try:
        new_cores_df = import_and_process_std_text_files(text_files_pattern, metadata_df)
    except Exception as e:
        print(f"Error processing text files: {e}")
        return
    
    if new_cores_df.empty:
        print("No cores processed successfully. Check file paths and formats.")
        return
    
    # Step 3: Save out standardized dataset
    else:
        try: 
            new_cores_df.to_csv(output_file_std, index=False)
            print(f"\n=== PROCESSING COMPLETE ===")
            print(f"Final dataset saved to: {output_file_std}")
            print(f"Cores with 550 kg/m³ horizon: {len(new_cores_df[~new_cores_df['depth_to_550'].isna()])}")
            print(f"Cores with 830 kg/m³ horizon: {len(new_cores_df[~new_cores_df['depth_to_830'].isna()])}")

        
        except Exception as e:
            print(f"Error creating final dataset: {e}")
            return

if __name__ == "__main__":
    process_std_main()