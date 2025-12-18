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
    
    # Remove extra index column
    df = df.drop(columns=['Unnamed: 0']) 
    
    #Fix elevation column
    df['Elevation'] = df['Elevation'].replace({'X': np.nan, 'NaN': np.nan, 'nan': np.nan})
    df["Elevation"] = df["Elevation"].astype(float)

    return df

def calculate_depth_to_density_std(depth, density, core_name, depth_dependence=True, do_plot=True):
    """Calculate depth where density first reaches target value using interpolation"""

    target_densities = [550, 830]
    colors = ['blue', 'red']
    depths_at_target_densities = [np.nan,np.nan]
    failed = [0,0]
    r2 = [np.nan, np.nan]

    models = def_models()

    # Remove NaN values
    valid_mask = ~(np.isnan(depth) | np.isnan(density))

    #If too few points, can't do fit
    if valid_mask.sum() < 10:
        
        failed = [1,1]
        print(" - Too few valid points ("+str(valid_mask.sum())+")")

        return depths_at_target_densities, failed, r2
    
    depth = np.array(depth[valid_mask])
    density = np.array(density[valid_mask])

    # Sort by depth
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
        ax.set_title(f'{core_name}')

    # target density out of range
    for idx, rho_i in enumerate(target_densities):

        # check if target density is within the range of measured densities and/or has more than 10 observatinos
        if rho_i < min_density or rho_i > max_density or len(depth_sorted) < 10:
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

        # Root finding to get depth at target density
        # If rho_i = 550, fits a curve between density.min and idx[density==600]
        # If rho_i = 830, fits a curve between idx[density==550]+1 and density.max
        # If depth_dependence fitting = false, fits curve to all data
        try: 
            
            if depth_dependence:
                if rho_i == 550:
                    fit_mask = density_sorted <= 600  # fit only to densities up to 600
                else:
                    fit_mask = density_sorted >= 550  # fit only to densities above 600
            else:
                fit_mask = np.ones_like(density_sorted, dtype=bool)  # use all data

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
            r2[idx] = round(1 - (ss_res / ss_tot),2)

            failed[idx] = 0            
    
        except Exception as e:
            print(f"Model fitting failed for {name} at {rho_i}: {e}")
    
        # try to find the depth at target density using the fit model
            
        # Root finding to get depth at target density    
        try:
            depth_at_rho = brentq(
                func_min,
                depth_sorted.min(),
                depth_sorted.max(),
                args=(*popt, rho_i)
            )
            
            depths_at_target_densities[idx] = round(depth_at_rho, 2)

        except Exception as e:
            print(f"  - Root finding failed for {rho_i}.")
            failed[idx] = 1
            pass
            

        if do_plot & (not failed[idx]):
            
            ax.plot(density_smooth, depth_smooth, '-', label=f'Model for {rho_i}', color=colors[idx])
            ax.plot(rho_i, depth_at_rho, 'x', color=colors[idx], markersize=10)

            ax.legend()
            ax.set_xlabel('Density (kg/m³)')
            ax.set_ylabel('Depth (m)')
            ax.set_xlim([250,950])
            
    
    if do_plot:
        plt.gca().invert_yaxis()
        plt.show()


    return depths_at_target_densities, failed, r2

def import_and_process_std_text_files(text_files, metadata_df, depth_dependence=True, do_plot=False):
    """Process all firn core std text files and calculate depth to 550 and 830"""
    
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
            depths_at_targets, failed, r2 = calculate_depth_to_density_std(data['depth'], data['density'], str(metadata_row['Name']), depth_dependence=depth_dependence, do_plot=do_plot)
            
            # Create core record
            core_record = {
                'core_name': str(metadata_row['Name']),
                'latitude': metadata_row['Latitude'],
                'longitude': metadata_row['Longitude'],
                'elevation': metadata_row.get('Elevation', np.nan),
                'year': metadata_row.get('Year', np.nan),
                'depth_to_550': depths_at_targets[0],
                'depth_to_830': depths_at_targets[1],
                'r2_550': r2[0],
                'r2_830': r2[1],
                'failed_fit_550': failed[0],
                'failed_fit_830': failed[1],
                'source': "PKM standardized cores",
                'method': metadata_row.get('Method', ''),
                'citation': metadata_row.get('Citation', ''),
                'measurement_count': len(data),
                'source_file': filename,
            }
            
            processed_cores.append(core_record)
            print(f"Processed {core_name}: {len(data)} measurements, depth_550={depths_at_targets[0]}, depth_830={depths_at_targets[1]}")
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    cores_df = pd.DataFrame(processed_cores)

    print(f"\nProcessed {len(cores_df)} cores successfully")
    
    return cores_df

def process_std_main(data_dir = "../data/PKM/", output_dir = "../data/PKM/processed/", domain='greenland', do_plot=False):
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
        text_files = glob.glob(text_files_pattern)
        print(f"Found {len(text_files)} text files")

        depth_dependence = True
        
        new_cores_df = import_and_process_std_text_files(text_files, metadata_df, depth_dependence, do_plot)
        new_cores_df["region"] = domain

    except Exception as e:
        print(f"Error processing text files: {e}")
        return
    
    # Step 3: Re-fit cores with all depths for poor 550 fits
    depth_dependence = False
    
    bad_fits = ['GGU163', 'Basin9', '6938', '7145', 'ACT_1']
    text_files_bad_fits = [f for f in text_files if any(bad in f for bad in bad_fits)]
    new_fits_for_bad_fits = import_and_process_std_text_files(text_files_bad_fits, metadata_df, depth_dependence=False, do_plot=False)

    for idx, row in new_fits_for_bad_fits.iterrows():

        if not np.isnan(row['depth_to_550']) or not np.isnan(row['depth_to_830']):
            replacement_row = new_fits_for_bad_fits[new_fits_for_bad_fits['core_name'] == row['core_name']]
            
            if not replacement_row.empty:
                print('replacing core: '+row['core_name'])
                for col in ['depth_to_550', 'depth_to_830', 'r2_550', 'r2_830', 'failed_fit_550', 'failed_fit_830']:
                    new_cores_df.loc[new_cores_df["core_name"]==row["core_name"], col] = replacement_row.iloc[0][col]

    # Step 4: Drop bad cores (cores where the depths are both np.nan, or cores where the data is too poor to trust)

    bad_cores = ['H1-1',  'H4-1', 'H2-1', 'GGU165', 'H3-1']

    clean_cores = new_cores_df[~new_cores_df['core_name'].isin(bad_cores)].copy(deep=True)
    clean_cores = clean_cores.reset_index(drop=True)

    drop_idxs = []
    for idx, row in clean_cores.iterrows():
        if np.isnan(row['depth_to_550']) and np.isnan(row['depth_to_830']):
            print('dropping core with no data: '+row['core_name'])
            drop_idxs.append(idx)
            #clean_cores = new_cores_df.drop(idx).copy(deep=True)
    clean_cores = clean_cores.drop(drop_idxs).reset_index(drop=True)

    
    if clean_cores.empty:
        print("No cores processed successfully. Check file paths and formats.")
        return
    
    # Step 5: Save out standardized dataset
    else:
        try: 
            clean_cores.to_csv(output_file, index=False)
            print(f"\n=== PROCESSING COMPLETE ===")
            print(f"Final dataset saved to: {output_file}")
            print(f"Cores with 550 kg/m³ horizon: {len(clean_cores[~clean_cores['depth_to_550'].isna()])}")
            print(f"Cores with 830 kg/m³ horizon: {len(clean_cores[~clean_cores['depth_to_830'].isna()])}")

        
        except Exception as e:
            print(f"Error creating final dataset: {e}")
            return

if __name__ == "__main__":
    process_std_main()