FC = mpifort
NetCDF = "netcdf4_DIR"
INCLUDE = -I/usr/local/apps/netcdf4/4.7.4/INTEL/2021.4/include
LIB2 = -L/usr/local/apps/netcdf4/4.7.4/INTEL/2021.4/lib -Wl,-rpath,/usr/local/apps/netcdf4/4.7.4/INTEL/2021.4/lib -lnetcdff -lnetcdf_c++ -lnetcdf
FFLAGS = -O3
EXTRA_0901 = -Wl,"--as-needed"

OBJECTS = \
  mainprogram_firnmodel.f90 \
  openNetCDF.f90 \
  ini_model.f90 \
  time_loop.f90 \
  output.f90 \
  subprogram.f90 

IMAU-FDM_np.x: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(FFLAGS) $(EXTRA_0901) $(LIB2) $(INCLUDE)

.f90.o:
	$(FC) -c $*.f90

