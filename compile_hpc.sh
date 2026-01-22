#!/bin/bash
# usage: compile_hpc.sh <profile>
# profile will be passed on to fpm's --profile flag. Default: release.

module reset

profile="$1"
if [[ -z "$profile" ]]; then
  profile="release"
fi

fpm=$(which fpm)

if [[ -z "$fpm" ]]; then
  module load conda
  # install fpm from conda
  conda create -c conda-forge --override-channels -n fpm -y fpm
  conda activate fpm
  fpm=$(which fpm)
  # create a symlink to fpm
  mkdir -p ~/.local/bin
  ln -s "$fpm" ~/.local/bin/fpm
  echo "Symlinked $fpm to ~/.local/bin/fpm"
  echo "Add ~/.local/bin/ to your PATH for convenience"
  conda deactivate
  module reset
fi

echo "Using fpm: $fpm"

module load prgenv/gnu &>/dev/null
module load hpcx-openmpi
module load netcdf4
module load nco
module load cdo
module load ncview
module load ncl
module load python3

NETCDF_LDFLAG=$(nc-config --libs)

(exec "$fpm" clean --all) || true
exec "$fpm" install --profile "$profile" --flag -ffree-line-length-512 --link-flag "$NETCDF_LDFLAG"
