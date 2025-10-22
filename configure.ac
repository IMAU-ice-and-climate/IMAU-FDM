AC_PREREQ([2.69])
AC_INIT([fdm], [1.0], [support@example.com])
AC_CONFIG_SRCDIR([source])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

# Check for Fortran compiler
AC_PROG_FC

# Set module directory variable
MOD_DIR="\$(top_builddir)/modules"
AC_SUBST([MOD_DIR])
#AC_CONFIG_COMMANDS([modules/mkdir], [true])
#AC_OUTPUT

# Check for DEBUG environment variable and set optimization flag
AC_MSG_CHECKING([for DEBUG mode])
if test "x$DEBUG" != "x"; then
  AC_MSG_RESULT([yes])
  DEBUG="true"
  OPTIMIZATION_FLAG="-Og"
  AC_MSG_NOTICE([Debug mode enabled, using optimization flag: $OPTIMIZATION_FLAG])
else
  AC_MSG_RESULT([no])
  DEBUG="false"
  OPTIMIZATION_FLAG="-O3"
  AC_MSG_NOTICE([Release mode, using optimization flag: $OPTIMIZATION_FLAG])
fi
AC_SUBST([OPTIMIZATION_FLAG])

# Check for NetCDF configuration tools and set NETCDF4_INCLUDE
AC_MSG_CHECKING([for NetCDF Fortran include flags])
if test "x$NETCDF4_INCLUDE" != "x"; then
  AC_MSG_RESULT([using environment variable: $NETCDF4_INCLUDE])
else
  # Try nc-config first
  if test "x`nc-config --fflags 2>/dev/null`" != "x"; then
    NETCDF4_INCLUDE=`nc-config --fflags`
    AC_MSG_RESULT([from nc-config: $NETCDF4_INCLUDE])
  # Fall back to nf-config
  elif test "x`nf-config --fflags 2>/dev/null`" != "x"; then
    NETCDF4_INCLUDE=`nf-config --fflags`
    AC_MSG_RESULT([from nf-config: $NETCDF4_INCLUDE])
  else
    AC_MSG_ERROR([Neither nc-config nor nf-config found. NetCDF support may be missing.])
  fi
fi
AC_SUBST([NETCDF4_INCLUDE])

# Check for NetCDF configuration tools and set NETCDF4_LIB
AC_MSG_CHECKING([for NetCDF Fortran library flags])
if test "x$NETCDF4_LIB" != "x"; then
  AC_MSG_RESULT([using environment variable: $NETCDF4_LIB])
else
  NETCDF4_LIB=""
  # Try nc-config for flibs first
  if test "x`nc-config --flibs 2>/dev/null`" != "x"; then
    NETCDF4_LIB=`nc-config --flibs`
    AC_MSG_NOTICE([Found nc-config --flibs: $NETCDF4_LIB])
  # Fall back to nf-config for flibs
  elif test "x`nf-config --flibs 2>/dev/null`" != "x"; then
    NETCDF4_LIB=`nf-config --flibs`
    AC_MSG_NOTICE([Found nf-config --flibs: $NETCDF4_LIB])
  else
    AC_MSG_ERROR([Neither nc-config --flibs nor nf-config --flibs found])
  fi
  
  # Add nc-config --libs to NETCDF4_LIB
  if test "x`nc-config --libs 2>/dev/null`" != "x"; then
    NC_LIBS=`nc-config --libs`
    NETCDF4_LIB="$NETCDF4_LIB $NC_LIBS"
    AC_MSG_NOTICE([Added nc-config --libs: $NC_LIBS])
  fi
  
  AC_MSG_RESULT([$NETCDF4_LIB])
fi
AC_SUBST([NETCDF4_LIB])

# Detect compiler and set flags
AC_MSG_CHECKING([Fortran compiler type])
case "$FC" in
  ifort*|*ifort)
    AC_MSG_RESULT([Intel Fortran (ifort)])
    FCFLAGS="-g -traceback -module \$(MOD_DIR)"
    FFLAGS="-Wall $OPTIMIZATION_FLAG -g -ffree-line-length-0 -fimplicit-none -fcheck=all -fbacktrace -I\$(MOD_DIR) $NETCDF4_INCLUDE $NETCDF4_LIB"
    AC_MSG_NOTICE([Intel Fortran compiler flags set])
    ;;
  gfortran*|*gfortran)
    AC_MSG_RESULT([GNU Fortran (gfortran)])
    FCFLAGS="-g -ffree-line-length-0 -fimplicit-none -fcheck=all -fbacktrace -J\$(MOD_DIR)"
    FFLAGS="$NETCDF4_INCLUDE $NETCDF4_LIB -I\$(MOD_DIR)"
    AC_MSG_NOTICE([GNU Fortran compiler flags set])
    ;;
  *)
    AC_MSG_RESULT([unknown: $FC])
    AC_MSG_WARN([Unknown Fortran compiler. Using default flags.])
    FCFLAGS="-g"
    FFLAGS="-g"
    ;;
esac

AC_MSG_CHECKING([for presence of mpifort])
AC_CHECK_PROG(MPIFORT, mpifort, yes)
AS_IF([test x"$MPIFORT" != x"yes"], [AC_MSG_ERROR([Please install mpifort before configuring.])])
# Use MPI Fortran compiler by default
ORIG_FC=$FC
FC=mpifort
AC_MSG_NOTICE([Setting mpifort as fortran compiler])
AC_SUBST([FF])

# Substitute variables
AC_SUBST([FCFLAGS])
AC_SUBST([FFLAGS])

# Output files
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_MSG_NOTICE([
========================================
Configuration Summary
========================================
Fortran compiler:     $FC ($ORIG_FC)
Module directory:     $MOD_DIR
Debug mode:           $DEBUG
Optimization flag:    $OPTIMIZATION_FLAG
NetCDF include flags: $NETCDF4_INCLUDE
NetCDF library flags: $NETCDF4_LIB

FCFLAGS: $FCFLAGS
FFLAGS:  $FFLAGS
========================================
])
