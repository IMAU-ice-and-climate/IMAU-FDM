# This script expects two environment variables to be set:
# NETCDF4_INCLUDE path to netcdf4-fortran include dir
# NETCDF4_LIB path to netcd4 shared libraries

# the DEBUG environment variable can optionally be set

# Define directory paths for source files, module files, and object files
SRC_DIR = source
MOD_DIR = modules
OBJ_DIR = objects

ifdef DEBUG
OPTIMIZATION_FLAG = -Og # optimize while keeping debug experience in mind
else
ifndef OPTIMIZATION_FLAG 
OPTIMIZATION_FLAG = -O3 # max optimization level
endif
endif

# Compiler and flags configuration
FC = mpifort # mpifort  # Specify the Fortran compiler

ifeq ($(OMPI_FC),ifort)
	FFLAGS = $(OPTIMIZATION_FLAG) -g -traceback -module $(MOD_DIR) $(NETCDF4_INCLUDE)  # Compilation flags: optimization, module directory, and NetCDF include
else #assumes gfortran 
	FFLAGS = -Wall $(OPTIMIZATION_FLAG) -g -ffree-line-length-0 -fimplicit-none -fcheck=all -fbacktrace -I $(MOD_DIR) $(NETCDF4_INCLUDE) # Compilation flags: optimization, module directory, and NetCDF include
endif

# List of source files
SRC_FILES = \
    model_settings.f90 \
    openNetCDF.f90 \
    output.f90 \
    grid_routines.f90 \
    water_physics.f90 \
    firn_physics.f90 \
    time_loop.f90 \
    initialise_variables.f90 \
    initialise_model.f90 \
    main.f90

# Generate a list of object files to be created from source files
OBJ_FILES = $(addprefix $(OBJ_DIR)/,$(patsubst %.f90,%.o,$(SRC_FILES)))

# Generate a list of corresponding module files
MOD_FILES = $(addprefix $(MOD_DIR)/,$(patsubst %.f90,%.mod,$(filter-out main.f90,$(SRC_FILES))))

# Define the target executable name
EXECUTABLE = imau-fdm.x

# Main target: compile the executable
all: $(EXECUTABLE)

# Rule to create the executable by linking object files
$(EXECUTABLE): $(OBJ_FILES)
	$(FC) -o $@ $^ $(FFLAGS) $(NETCDF4_LIB)

# Rule to compile each source file into an object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# Create the object directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Rule to generate empty module files (used to resolve dependencies)
$(MOD_DIR)/%.mod: $(SRC_DIR)/%.f90 | $(MOD_DIR)
	$(FC) $(FFLAGS) -c -o /dev/null $<

# Create the module directory if it doesn't exist
$(MOD_DIR):
	mkdir -p $(MOD_DIR)

# Clean target: remove object and module files, and the executable
clean:
	rm -f $(OBJ_DIR)/* $(MOD_DIR)/* $(EXECUTABLE)

# Declare "all" and "clean" as phony targets (no actual files associated)
.PHONY: all clean

