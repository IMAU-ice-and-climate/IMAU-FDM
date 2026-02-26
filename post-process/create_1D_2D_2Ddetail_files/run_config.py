"""
Auto-detecting run configuration for IMAU-FDM post-processing.

This module automatically extracts configuration from:
- NetCDF output files (domain, time range, timestep, variables)
- Model settings files (spinup, grid dimensions, output frequencies)
- Reference files (mask, pointlist)

Usage:
    from run_config import RunConfig

    config = RunConfig(
        output_dir='/path/to/output/',
        reference_dir='/path/to/reference/'
    )
    # All parameters are now auto-detected
    print(config.domain)
    print(config.grid_shape)

    Created by Elizabeth Case and Claude Code in January 2026
"""

import re
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple
import numpy as np


@dataclass
class RunConfig:
    """
    Auto-detecting configuration for IMAU-FDM post-processing.

    Parameters
    ----------
    output_dir : str or Path
        Directory containing model output files (1D, 2D, 2Ddetail)
    reference_dir : str or Path, optional
        Directory containing reference files (masks, pointlists).
        If None, attempts to find based on domain.
    ms_dir : str or Path, optional
        Directory containing model_settings files.
        If None, looks for 'ms_files' sibling to output_dir.
    processed_output_dir : str or Path, optional
        Directory for processed output files.
        If None, creates 'post-process' sibling to output_dir.

    Attributes (auto-detected)
    --------------------------
    domain : str
        Domain name (e.g., 'FGRN055', 'ANT27')
    model_start : datetime
        Model start date
    model_end : datetime
        Model end date
    timestep_1d : int
        Timestep for 1D files in seconds
    timestep_2d : int
        Timestep for 2D profile files in seconds
    timestep_2ddetail : int
        Timestep for 2Ddetail files in seconds
    spinup_years : int
        Number of spinup years
    total_years : int
        Total simulation years
    grid_shape : tuple
        Grid dimensions (nlat, nlon)
    n_points : int
        Number of grid points with data
    variables_1d : list
        Available variables in 1D files
    variables_2d : list
        Available variables in 2D files
    variables_2ddetail : list
        Available variables in 2Ddetail files
    detail_layers : int
        Number of layers in 2Ddetail files
    detail_thickness : float
        Layer thickness in 2Ddetail files (m)
    """
    output_dir: Path
    reference_dir: Optional[Path] = None
    ms_dir: Optional[Path] = None
    processed_output_dir: Optional[Path] = None

    # Auto-detected attributes
    domain: str = field(default=None, init=False)
    model_start: datetime = field(default=None, init=False)
    model_end: datetime = field(default=None, init=False)
    timestep_1d: int = field(default=None, init=False)
    timestep_2d: int = field(default=None, init=False)
    timestep_2ddetail: int = field(default=None, init=False)
    spinup_years: int = field(default=None, init=False)
    total_years: int = field(default=None, init=False)
    grid_shape: Tuple[int, int] = field(default=None, init=False)
    n_points: int = field(default=None, init=False)
    n_timesteps_1d: int = field(default=None, init=False)
    n_timesteps_2d: int = field(default=None, init=False)
    n_timesteps_2ddetail: int = field(default=None, init=False)
    variables_1d: List[str] = field(default_factory=list, init=False)
    variables_2d: List[str] = field(default_factory=list, init=False)
    variables_2ddetail: List[str] = field(default_factory=list, init=False)
    detail_layers: int = field(default=None, init=False)
    detail_thickness: float = field(default=None, init=False)
    profile_layers: int = field(default=None, init=False)

    # File patterns (set after domain detection)
    _file_pattern_1d: str = field(default=None, init=False)
    _file_pattern_2d: str = field(default=None, init=False)
    _file_pattern_2ddetail: str = field(default=None, init=False)

    def __post_init__(self):
        """Auto-detect all configuration after initialization."""
        self.output_dir = Path(self.output_dir)

        if not self.output_dir.exists():
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")

        # Set up directories
        self._setup_directories()

        # Auto-detect from NetCDF files
        self._detect_from_netcdf()

        # Auto-detect from model settings
        self._detect_from_model_settings()

        # Detect grid from mask as fallback
        self._detect_grid_from_mask()

        # Set file patterns
        self._set_file_patterns()

        # Load pointlist to get n_points
        self._detect_from_pointlist()

    def _setup_directories(self):
        """Set up reference, model settings, and output directories."""
        # Reference directory
        if self.reference_dir is None:
            # Try common locations
            possible_refs = [
                Path('/home/nld4814/perm/code/IMAU-FDM/reference'),
                self.output_dir.parent.parent / 'reference',
                self.output_dir.parent / 'reference',
            ]
            for ref in possible_refs:
                if ref.exists():
                    self.reference_dir = ref
                    break
        else:
            self.reference_dir = Path(self.reference_dir)

        # Model settings directory
        if self.ms_dir is None:
            possible_ms = [
                self.output_dir.parent / 'ms_files',
                self.output_dir / 'ms_files',
            ]
            for ms in possible_ms:
                if ms.exists():
                    self.ms_dir = ms
                    break
        else:
            self.ms_dir = Path(self.ms_dir)

        # Processed output directory
        if self.processed_output_dir is None:
            self.processed_output_dir = self.output_dir.parent / 'post-process'
        else:
            self.processed_output_dir = Path(self.processed_output_dir)

    def _detect_from_netcdf(self):
        """Extract configuration from NetCDF output files."""
        import xarray as xr

        # Find a sample 1D file
        sample_1d = self._find_sample_file('1D')
        if sample_1d:
            with xr.open_dataset(sample_1d) as ds:
                # Domain
                self.domain = ds.attrs.get('domain', 'UNKNOWN')

                # Time range
                start_str = ds.attrs.get('model_start_datetime', '')
                end_str = ds.attrs.get('model_end_datetime', '')
                if start_str:
                    self.model_start = datetime.fromisoformat(start_str.replace('T', ' ').split('.')[0])
                if end_str:
                    self.model_end = datetime.fromisoformat(end_str.replace('T', ' ').split('.')[0])

                # Timestep
                ts_str = ds.attrs.get('timestep_length_of_dimension_in_seconds', '86400')
                self.timestep_1d = int(ts_str)

                # Number of timesteps
                time_dim = [d for d in ds.sizes if 'ind_t' in d or 'time' in d.lower()]
                if time_dim:
                    self.n_timesteps_1d = ds.sizes[time_dim[0]]

                # Variables (exclude coordinates)
                self.variables_1d = [v for v in ds.data_vars
                                     if v not in ['lat', 'lon', 'time', 'rlat', 'rlon']]

        # Find a sample 2D file
        sample_2d = self._find_sample_file('2D')
        if sample_2d:
            with xr.open_dataset(sample_2d) as ds:
                ts_str = ds.attrs.get('timestep_length_of_dimension_in_seconds', '2592000')
                self.timestep_2d = int(ts_str)

                # Get dimensions
                if 'layer' in ds.sizes:
                    self.profile_layers = ds.sizes['layer']
                time_dim = [d for d in ds.sizes if 'ind_t' in d or 'time' in d.lower()]
                if time_dim:
                    self.n_timesteps_2d = ds.sizes[time_dim[0]]

                self.variables_2d = [v for v in ds.data_vars]

        # Find a sample 2Ddetail file
        sample_2ddetail = self._find_sample_file('2Ddetail')
        if sample_2ddetail:
            with xr.open_dataset(sample_2ddetail) as ds:
                ts_str = ds.attrs.get('timestep_length_of_dimension_in_seconds', '864000')
                self.timestep_2ddetail = int(ts_str)

                # Get dimensions
                if 'layer' in ds.sizes:
                    self.detail_layers = ds.sizes['layer']
                time_dim = [d for d in ds.sizes if 'ind_t' in d or 'time' in d.lower()]
                if time_dim:
                    self.n_timesteps_2ddetail = ds.sizes[time_dim[0]]

                # Get layer thickness from depth or dz variable
                if 'dz' in ds:
                    self.detail_thickness = float(ds['dz'].values[0])
                elif 'depth' in ds and ds['depth'].ndim == 1:
                    depths = ds['depth'].values
                    if len(depths) > 1:
                        self.detail_thickness = float(depths[1] - depths[0])

                self.variables_2ddetail = [v for v in ds.data_vars]

    def _find_sample_file(self, file_type: str) -> Optional[Path]:
        """Find a sample file of the given type."""
        patterns = {
            '1D': '*_1D_*.nc',
            '2D': '*_2D_[0-9]*.nc',  # Exclude 2Ddetail
            '2Ddetail': '*_2Ddetail_*.nc',
        }

        pattern = patterns.get(file_type, '*.nc')
        files = list(self.output_dir.glob(pattern))

        # For 2D, exclude 2Ddetail files
        if file_type == '2D':
            files = [f for f in files if '2Ddetail' not in f.name]

        return files[0] if files else None

    def _detect_from_model_settings(self):
        """Extract configuration from model_settings file."""
        if self.ms_dir is None or not self.ms_dir.exists():
            return

        # Find a model settings file
        ms_files = list(self.ms_dir.glob('model_settings_*.txt'))
        if not ms_files:
            return

        ms_file = ms_files[0]

        # Parse the file - format is "value ! key; description" or "value ! key, description"
        settings = {}
        with open(ms_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('!-'):
                    continue

                # Parse "value ! comment" format
                if '!' in line:
                    value_part, comment = line.split('!', 1)
                    value = value_part.strip()

                    # Extract key from comment (handles both "key;" and "key," formats)
                    comment = comment.strip()
                    key_match = re.match(r'(\w+)[;,]', comment)
                    if key_match:
                        key = key_match.group(1)
                        settings[key] = value

        # Look for spinup years (second occurrence of nyears)
        lines = open(ms_file).readlines()
        nyears_count = 0
        for line in lines:
            if '!' in line and 'nyears' in line:
                value = line.split('!')[0].strip()
                try:
                    val = int(value)
                    nyears_count += 1
                    if nyears_count == 1:
                        self.total_years = val
                    elif nyears_count == 2:
                        self.spinup_years = val
                except ValueError:
                    pass

        # Grid dimensions from settings
        if 'numLats' in settings and 'numLons' in settings:
            try:
                nlat = int(settings['numLats'])
                nlon = int(settings['numLons'])
                self.grid_shape = (nlat, nlon)
            except ValueError:
                pass

        # Detail settings
        if 'detlayers' in settings and self.detail_layers is None:
            try:
                self.detail_layers = int(settings['detlayers'])
            except ValueError:
                pass
        if 'detthick' in settings and self.detail_thickness is None:
            try:
                self.detail_thickness = float(settings['detthick'])
            except ValueError:
                pass

    def _detect_grid_from_mask(self):
        """Detect grid shape from mask file if not already set."""
        if self.grid_shape is not None:
            return

        mask_path = self.get_mask_path()
        if mask_path and mask_path.exists():
            import xarray as xr
            with xr.open_dataset(mask_path) as ds:
                # Get shape from lat/lon arrays
                if 'lat' in ds:
                    shape = ds['lat'].shape
                    if len(shape) == 2:
                        self.grid_shape = shape

    def _detect_from_pointlist(self):
        """Get number of points from pointlist."""
        pointlist = self.get_pointlist_path()
        if pointlist and pointlist.exists():
            with open(pointlist, 'r') as f:
                self.n_points = sum(1 for line in f if line.strip())

    def _set_file_patterns(self):
        """Set file name patterns based on detected domain."""
        if self.domain:
            # Detect the actual pattern from existing files
            sample = self._find_sample_file('1D')
            if sample:
                # Extract pattern: e.g., "FGRN055_era055_1D_{}.nc"
                name = sample.name
                # Replace the point number with {}
                pattern = re.sub(r'_1D_\d+\.nc$', '_1D_{}.nc', name)
                self._file_pattern_1d = pattern
                self._file_pattern_2d = pattern.replace('_1D_', '_2D_')
                self._file_pattern_2ddetail = pattern.replace('_1D_', '_2Ddetail_')

    # Path getters
    def get_mask_path(self) -> Optional[Path]:
        """Get path to mask file for this domain."""
        if self.reference_dir and self.domain:
            mask_path = self.reference_dir / self.domain / f'{self.domain}_Masks.nc'
            if mask_path.exists():
                return mask_path
        return None

    def get_pointlist_path(self) -> Optional[Path]:
        """Get path to pointlist file for this domain."""
        if self.reference_dir and self.domain:
            domain_dir = self.reference_dir / self.domain
            if domain_dir.exists():
                # Try common patterns
                patterns = [
                    f'IN_ll_{self.domain}.txt',
                    f'IN_ll_{self.domain}_*.txt',
                ]
                for pattern in patterns:
                    files = list(domain_dir.glob(pattern))
                    if files:
                        return files[0]
        return None

    def get_1d_file(self, point_num: int) -> Path:
        """Get path to 1D file for a given point number."""
        if self._file_pattern_1d:
            return self.output_dir / self._file_pattern_1d.format(point_num)
        return self.output_dir / f'{self.domain}_1D_{point_num}.nc'

    def get_2d_file(self, point_num: int) -> Path:
        """Get path to 2D file for a given point number."""
        if self._file_pattern_2d:
            return self.output_dir / self._file_pattern_2d.format(point_num)
        return self.output_dir / f'{self.domain}_2D_{point_num}.nc'

    def get_2ddetail_file(self, point_num: int) -> Path:
        """Get path to 2Ddetail file for a given point number."""
        if self._file_pattern_2ddetail:
            return self.output_dir / self._file_pattern_2ddetail.format(point_num)
        return self.output_dir / f'{self.domain}_2Ddetail_{point_num}.nc'

    def get_output_filename(self, var_name: str, file_type: str = '1D') -> str:
        """Generate output filename for processed data."""
        start_year = self.model_start.year if self.model_start else 'XXXX'
        end_year = self.model_end.year if self.model_end else 'XXXX'
        return f'FDM_{var_name}_{self.domain}_{start_year}-{end_year}_{file_type}.nc'

    def get_fractional_year(self, timestep_idx: int, file_type: str = '1D') -> float:
        """
        Convert timestep index to fractional year.

        Parameters
        ----------
        timestep_idx : int
            Zero-based timestep index
        file_type : str
            File type ('1D', '2D', '2Ddetail')

        Returns
        -------
        float
            Fractional year
        """
        if self.model_start is None:
            return float(timestep_idx)

        # Get timestep in seconds
        if file_type == '1D':
            dt = self.timestep_1d or 86400
        elif file_type == '2D':
            dt = self.timestep_2d or 2592000
        else:  # 2Ddetail
            dt = self.timestep_2ddetail or 864000

        # Calculate fractional year
        # Start from model start, add timestep_idx * dt seconds
        start_year = self.model_start.year
        start_doy = self.model_start.timetuple().tm_yday

        # Starting fractional year
        days_in_year = 366 if self.model_start.year % 4 == 0 else 365
        start_frac = start_year + (start_doy - 1) / days_in_year

        # Add timestep contribution (dt in days / days per year)
        dt_days = dt / 86400.0
        return start_frac + (timestep_idx * dt_days) / 365.25

    def get_spinup_end_index(self, file_type: str = '1D') -> int:
        """
        Get the timestep index where spinup ends.

        Returns
        -------
        int
            Index of first post-spinup timestep
        """
        if self.spinup_years is None or self.model_start is None:
            return 0

        # Get timestep in days
        if file_type == '1D':
            dt_days = (self.timestep_1d or 86400) / 86400.0
        elif file_type == '2D':
            dt_days = (self.timestep_2d or 2592000) / 86400.0
        else:
            dt_days = (self.timestep_2ddetail or 864000) / 86400.0

        # Spinup covers spinup_years
        spinup_days = self.spinup_years * 365.25
        return int(spinup_days / dt_days)

    def summary(self) -> str:
        """Return a summary of the detected configuration."""
        lines = [
            "=" * 60,
            "IMAU-FDM Run Configuration",
            "=" * 60,
            f"Domain:              {self.domain}",
            f"Model period:        {self.model_start} to {self.model_end}",
            f"Total years:         {self.total_years}",
            f"Spinup years:        {self.spinup_years}",
            f"Grid shape:          {self.grid_shape}",
            f"Number of points:    {self.n_points}",
            "",
            "File types:",
            f"  1D:       {self.n_timesteps_1d} timesteps @ {self.timestep_1d}s",
            f"  2D:       {self.n_timesteps_2d} timesteps @ {self.timestep_2d}s ({self.profile_layers} layers)",
            f"  2Ddetail: {self.n_timesteps_2ddetail} timesteps @ {self.timestep_2ddetail}s ({self.detail_layers} layers @ {self.detail_thickness}m)",
            "",
            "Directories:",
            f"  Output:     {self.output_dir}",
            f"  Reference:  {self.reference_dir}",
            f"  Settings:   {self.ms_dir}",
            f"  Processed:  {self.processed_output_dir}",
            "",
            "Reference files:",
            f"  Mask:       {self.get_mask_path()}",
            f"  Pointlist:  {self.get_pointlist_path()}",
            "",
            f"Variables (1D):       {', '.join(self.variables_1d[:5])}{'...' if len(self.variables_1d) > 5 else ''}",
            f"Variables (2D):       {', '.join(self.variables_2d)}",
            f"Variables (2Ddetail): {', '.join(self.variables_2ddetail)}",
            "=" * 60,
        ]
        return '\n'.join(lines)

    def __repr__(self):
        return f"RunConfig(domain='{self.domain}', grid={self.grid_shape}, points={self.n_points})"


def load_pointlist(pointlist_path: Path) -> np.ndarray:
    """
    Load pointlist file and return as structured array.

    Parameters
    ----------
    pointlist_path : Path
        Path to pointlist file

    Returns
    -------
    np.ndarray
        Structured array with columns: lon, lat, ..., rlat_idx, rlon_idx
    """
    data = np.loadtxt(pointlist_path, delimiter=',')
    return data


def load_mask(mask_path: Path):
    """
    Load mask file and return coordinates and mask.

    Parameters
    ----------
    mask_path : Path
        Path to mask NetCDF file

    Returns
    -------
    dict
        Dictionary with 'lat', 'lon', 'mask', 'rlat', 'rlon'
    """
    import xarray as xr

    with xr.open_dataset(mask_path) as ds:
        result = {
            'lat': ds['lat'].values,
            'lon': ds['lon'].values,
            'rlat': ds['rlat'].values if 'rlat' in ds else np.arange(ds.dims.get('rlat', ds.dims.get('y', 0))),
            'rlon': ds['rlon'].values if 'rlon' in ds else np.arange(ds.dims.get('rlon', ds.dims.get('x', 0))),
        }

        # Find mask variable
        for mask_name in ['IceMask', 'icemask', 'mask', 'Icemask']:
            if mask_name in ds:
                result['mask'] = ds[mask_name].values
                break

        return result


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        output_dir = sys.argv[1]
    else:
        output_dir = '/home/nld4814/scratch/run_FGRN055-era055_1939-2023/output/'

    config = RunConfig(output_dir=output_dir)
    print(config.summary())
