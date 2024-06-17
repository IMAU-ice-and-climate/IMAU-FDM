#!/bin/ksh -x

## slices decadal racmo output into year files                  ##
## ------------------------------------------------------------ ##

project_name="FGRN055-era055"
raw_dir="$HPCPERM/${project_name}/raw/historical/"
fname_extra=".FGRN055.BN_RACMO2.3p2_FGRN055.3H.nc"

## will be one of: precip snowfall evap tskin sndiv snowmelt    ##
## ------------------------------------------------------------ ##
varlist=$1

## sets number of data points in each year (#d/yr * #pts/d)     ##
## TK: update to 2023?
## ------------------------------------------------------------ ##
set -A yearlist 1957 1961 1971 1981 1991 2001 2011
set -A nyears   4    10   10   10   10   10   10

set -A Y57list 736 2920 2920 2928
set -A Y61list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
set -A Y71list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
set -A Y81list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
set -A Y91list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928 
set -A Y01list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
set -A Y11list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928