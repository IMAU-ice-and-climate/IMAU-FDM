SRC_DIR = source

bin_PROGRAMS = fdm
fdm_SOURCES = \
    $(SRC_DIR)/model_settings.f90 \
    $(SRC_DIR)/openNetCDF.f90 \
    $(SRC_DIR)/output.f90 \
    $(SRC_DIR)/grid_routines.f90 \
    $(SRC_DIR)/water_physics.f90 \
    $(SRC_DIR)/firn_physics.f90 \
    $(SRC_DIR)/time_loop.f90 \
    $(SRC_DIR)/initialise_variables.f90 \
    $(SRC_DIR)/initialise_model.f90 \
    $(SRC_DIR)/main.f90

# The final flags passed to the compiler will be "$(AM_FCFLAGS) $(FCFLAGS)"
# So we can pass additional flags by overriding $(AM_FCFLAGS)
AM_FCFLAGS = $(FFLAGS)

all-local: $(MOD_DIR)

$(MOD_DIR):
	$(MKDIR_P) $(MOD_DIR)

# Make all object files depend on the module directory existing
$(fdm_OBJECTS): | $(MOD_DIR)

# Clean module files and directory
clean-local:
	rm -rf $(MOD_DIR)

CLEANFILES = $(MOD_DIR)/*.mod
