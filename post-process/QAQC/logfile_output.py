import re
import numpy as np
import matplotlib.pyplot as plt

def extract_spinup_data(logfile_path):
    """
    Extracts spinup data (variable data)
    
    Returns a dictionary with spinup number as key and variable dict as value.
    """
    spinup_data = {}
    
    with open(logfile_path, 'r') as f:
        content = f.read()
    
    # Split by "After spin-up #" sections
    sections = re.split(r'After spin-up #\s+(\d+)', content)
    
    # Process each spinup section (skip first element which is pre-spinup content)
    for i in range(1, len(sections), 2):
        spinup_num = int(sections[i])
        section_text = sections[i+1]
        
        # Extract variables using regex
        data = {}
        
        # Extract Rho(200)
        rho_match = re.search(r'Rho\(200\)\s*=\s*([\d.E+-]+)', section_text)
        if rho_match:
            data['Rho_200'] = float(rho_match.group(1))
        
        # Extract T(200)
        t_match = re.search(r'T\(200\)\s*=\s*([\d.E+-]+)', section_text)
        if t_match:
            data['T_200'] = float(t_match.group(1))
        
        # Extract Year(200)
        year_match = re.search(r'Year\(200\)\s*=\s*([-\d.E+-]+)', section_text)
        if year_match:
            data['Year_200'] = float(year_match.group(1))
        
        # Extract Ind_z_surf
        ind_match = re.search(r'Ind_z_surf\s*=\s*([\d]+)', section_text)
        if ind_match:
            data['Ind_z_surf'] = int(ind_match.group(1))
        
        # Extract h_surf
        hsurf_match = re.search(r'h_surf\s*=\s*([-\d.E+-]+)', section_text)
        if hsurf_match:
            data['h_surf'] = float(hsurf_match.group(1))
        
        # Extract FAC
        fac_match = re.search(r'FAC\s*=\s*([\d.E+-]+)', section_text)
        if fac_match:
            data['FAC'] = float(fac_match.group(1))
        
        # Extract IceMass
        icemass_match = re.search(r'IceMass\s*=\s*([\d.E+-]+)', section_text)
        if icemass_match:
            data['IceMass'] = float(icemass_match.group(1))
        
        spinup_data[spinup_num] = data
    
    return spinup_data

def plot_spinup_variables(spinup_data, point_num, figsize=(15, 10)):
    """
    Plot all spinup variables vs spinup number.
    """
    # Sort by spinup number
    spinup_nums = sorted(spinup_data.keys())
    
    # Extract variable names (from first entry)
    var_names = list(spinup_data[spinup_nums[0]].keys())
    
    # Create subplots
    n_vars = len(var_names)
    n_cols = 3
    n_rows = int(np.ceil(n_vars / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()
    
    # Plot each variable
    for idx, var_name in enumerate(var_names):
        values = [spinup_data[num].get(var_name, np.nan) for num in spinup_nums]
        
        axes[idx].plot(spinup_nums, values, 'o-', linewidth=2, markersize=4)
        axes[idx].set_xlabel('Spin-up #', fontsize=10)
        axes[idx].set_ylabel(var_name, fontsize=10)
        axes[idx].set_title(var_name, fontsize=11, fontweight='bold')
        axes[idx].grid(True, alpha=0.3)
    
    # Hide extra subplots if any
    for idx in range(n_vars, len(axes)):
        axes[idx].set_visible(False)
    
    plt.title(point_num)
    plt.tight_layout()
    return fig

# Plot individual variables with more detail
def plot_individual_variable(spinup_data, var_name, point_num):
    spinup_nums = sorted(spinup_data.keys())
    values = [spinup_data[num].get(var_name, np.nan) for num in spinup_nums]
    
    plt.figure(figsize=(10, 6))
    plt.plot(spinup_nums, values, 'o-', linewidth=2, markersize=6)
    plt.xlabel('Spin-up #', fontsize=12)
    plt.ylabel(var_name, fontsize=12)
    plt.title(f'{var_name} for {point_num}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    return plt.gcf()

