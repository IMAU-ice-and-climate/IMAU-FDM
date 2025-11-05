import pandas as pd
import numpy as np
import os
import glob
import warnings
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, brentq
warnings.filterwarnings('ignore')
from process_sumup_cores import def_models

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
    
    return df

def calculate_depth_to_density(depth, density, target_density, do_plot=False):
    """Calculate depth where density first reaches target value using interpolation"""
    #TODO: update this interpolation function to match the sumup one

    models = def_models()

    failed=0
    depth_at_rho = np.nan
    best_r2 = -np.inf
    best_model = None

    # Remove NaN values
    valid_mask = ~(np.isnan(depth) | np.isnan(density))

    if valid_mask.sum() < 2:
        depth_at_rho = np.nan
        failed = 1
        return depth_at_rho, failed, best_r2
    
    depth = np.array(depth)
    density = np.array(density)
    
    # Sort by depth
    sort_idx = np.argsort(depth)
    depth_sorted = depth[sort_idx]
    density_sorted = density[sort_idx]
    
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
    # # Find first point where target density is reached or exceeded
    # idx = np.where(densities_sorted >= target_density)[0]
    # if len(idx) == 0:
    #     return np.nan
    
    # first_idx = idx[0]
    
    # # If first measurement already exceeds target
    # if first_idx == 0:
    #     return depths_sorted[0]
    
    # # Linear interpolation between the two points
    # prev_depth = depths_sorted[first_idx - 1]
    # prev_density = densities_sorted[first_idx - 1]
    # curr_depth = depths_sorted[first_idx]
    # curr_density = densities_sorted[first_idx]
    
    # # Interpolate
    # if curr_density == prev_density:
    #     return prev_depth
    
    # interpolated_depth = prev_depth + (target_density - prev_density) * (curr_depth - prev_depth) / (curr_density - prev_density)
    
    # return round(interpolated_depth, 2)

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
            depth_550, failed_550, r2_550 = calculate_depth_to_density(data['depth'], data['density'], 550)
            depth_830, failed_830, r2_830 = calculate_depth_to_density(data['depth'], data['density'], 830)
            
            # Create core record
            core_record = {
                'core_name': str(metadata_row['Name']),
                'latitude': metadata_row['Latitude'],
                'longitude': metadata_row['Longitude'],
                'elevation': metadata_row.get('Elevation', np.nan),
                'year': metadata_row.get('Year', np.nan),
                'depth_to_550': depth_550,
                'depth_to_830': depth_830,
                'r2_550': r2_550,
                'r2_830': r2_830,
                'failed_fit_550': failed_550,
                'failed_fit_830': failed_830,
                'source': "PKM standardized cores",
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

def process_std_main(data_dir = "../data/PKM/", output_dir = "../data/PKM/processed/", domain='greenland'):
    """Main processing function"""

    metadata_file = data_dir+"cleaned-std-metadata.csv"
    text_files_pattern = data_dir+"source/*.txt"  # or "*.std.txt" or "/path/to/cores/*.txt"
    output_file = f"{data_dir}processed/PKM_550_830_density_depths_{domain}.csv"
    
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
            new_cores_df.to_csv(output_file, index=False)
            print(f"\n=== PROCESSING COMPLETE ===")
            print(f"Final dataset saved to: {output_file}")
            print(f"Cores with 550 kg/m³ horizon: {len(new_cores_df[~new_cores_df['depth_to_550'].isna()])}")
            print(f"Cores with 830 kg/m³ horizon: {len(new_cores_df[~new_cores_df['depth_to_830'].isna()])}")

        
        except Exception as e:
            print(f"Error creating final dataset: {e}")
            return

if __name__ == "__main__":
    process_std_main()