# To build an image from this file using docker or podman, from this directory:
# docker build --platform=linux/amd64 -t imau-fdm -f Containerfile .
FROM ubuntu:jammy
WORKDIR /imau-fdm
COPY . .
# Install requirements
RUN apt update -y && DEBIAN_FRONTEND=noninteractive && \
    apt install -y --no-install-recommends make mpich gfortran libmpich-dev \
    ca-certificates curl jq m4 \
    hdf5-tools hdf5-helpers libhdf5-dev libhdf5-serial-dev libnetcdf-dev # netcdf4-fortran dependencies
# Download netcdf-fortran from github, then build and install it
RUN export NETCDF_URL=$(curl https://api.github.com/repos/Unidata/netcdf-fortran/releases/latest | jq -r '.tarball_url') && \
    curl -L "$NETCDF_URL" --output netcdf.tar.gz && \
    tar -xzf netcdf.tar.gz && \
    cd Unidata-netcdf* && \
    ./configure && make && make install && ldconfig && \
    cd .. && rm -rf Unidata-netcdf* && \
    nf-config --all
# Make imau-fdm
RUN DEBUG=true NETCDF4_LIB=$(nc-config --flibs) NETCDF4_INCLUDE=$(nc-config --fflags) make && rm -rf objects/** modules/**
CMD ./imau-fdm.x
