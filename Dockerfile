# THIS FILE CONTAINS THE NECESSARY ELEMENTS TO RUN MEDSLIK-II in a container

# syntax=docker/dockerfile:1

ARG PYTHON_VERSION=3.11.4
FROM python:${PYTHON_VERSION}-slim as base

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

#Set up the Home dor as the Medslik Directory
WORKDIR /Medslik-II

# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/Medslik-II" \
    --uid "${UID}" \
    appuser

#Allow appuser to navigate and read scripts
RUN chown -R appuser:appuser /Medslik-II

# Download dependencies as a separate step to take advantage of Docker's caching.
# Leverage a cache mount to /root/.cache/pip to speed up subsequent builds.
# Leverage a bind mount to requirements.txt to avoid having to copy them into
# into this layer.
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    python -m pip install -r requirements.txt

RUN apt-get update \
  && apt-get install -yq --no-install-recommends \
  build-essential \
  curl \
  fuse \
  gfortran \
  g++ \
  git \
  gnupg \
  gnupg2 \
  keychain \
  libcurl4-openssl-dev \
  libfuse-dev \
  liblapack-dev \
  libssl-dev \
  locate \
  lsb-release \
  make \
  m4 \
  nano \
  rsync \
  tzdata \
  tini \
  unzip \
  vim \
  wget \
  zip
  
# build netcdf with gcc and g-fortran
ENV CC=gcc
ENV FC=gfortran

# set library location
ENV PREFIXDIR=/usr/local

WORKDIR /

## get zlib
ENV ZLIB_VERSION=zlib-1.3.1
RUN wget https://zlib.net/${ZLIB_VERSION}.tar.gz && tar -xvzf ${ZLIB_VERSION}.tar.gz
RUN cd ${ZLIB_VERSION} \
    && ./configure --prefix=${PREFIXDIR} \
    && make install
WORKDIR /
RUN rm -rf ${ZLIB_VERSION}.tar.gz ${ZLIB_VERSION}

## get hdf5-1.8
ENV HDF518_VERSION=hdf5-1.8.21
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/${HDF518_VERSION}/src/${HDF518_VERSION}.tar.gz && tar -xvzf ${HDF518_VERSION}.tar.gz
RUN cd ${HDF518_VERSION} \
    && ./configure --with-zlib=${PREFIXDIR} --prefix=${PREFIXDIR} --enable-hl --enable-shared\
    && make \
    && make install
WORKDIR /
RUN rm -rf ${HDF518_VERSION}.tar.gz ${HDF518_VERSION}

## get netcdf-c
ENV NETCDFC_VERSION=4.8.0
RUN wget https://github.com/Unidata/netcdf-c/archive/v${NETCDFC_VERSION}.tar.gz && tar -xvzf v${NETCDFC_VERSION}.tar.gz
RUN cd netcdf-c-${NETCDFC_VERSION} \
    && CPPFLAGS=-I${PREFIXDIR}/include LDFLAGS=-L${PREFIXDIR}/lib ./configure --prefix=${PREFIXDIR} --enable-netcdf-4 --enable-shared --enable-dap \
    && make install
WORKDIR /
RUN rm -rf v${NETCDFC_VERSION}.tar.gz netcdf-fortran-${NETCDFC_VERSION}

# Set these flags because some NETCDF libraries do not allow some numerical warnings. As they are not error, these f;ags do not interfer with the results
ENV FCFLAGS="-w -fallow-argument-mismatch -O2"
ENV FFLAGS="-w -fallow-argument-mismatch -O2"

## get netcdf-fortran
ENV NETCDFFORTRAN_VERSION=4.6.0
RUN wget https://github.com/Unidata/netcdf-fortran/archive/v${NETCDFFORTRAN_VERSION}.tar.gz && tar -xvzf v${NETCDFFORTRAN_VERSION}.tar.gz
RUN cd netcdf-fortran-${NETCDFFORTRAN_VERSION} \
    && CPPFLAGS=-I${PREFIXDIR}/include LDFLAGS=-L${PREFIXDIR}/lib ./configure --prefix=${PREFIXDIR} \
    # && make check \
    && make install
WORKDIR /
RUN rm -rf v${NETCDFFORTRAN_VERSION}.tar.gz netcdf-fortran-${NETCDFFORTRAN_VERSION}

#Create a env variable to indicate where the netcdf libaries are saved. Otherwise, executables are not able to run
ENV LD_LIBRARY_PATH=/usr/local/lib

# Switch to the non-privileged user to run the application.
USER appuser

#Change to home directory again
WORKDIR /Medslik-II

# Copy the source code into the container.
COPY --chown=appuser:appuser . .

# Expose the port that the application listens on.
EXPOSE 8000

# Run the application.
#This command allows to navigate and test the code in the docker machine. This will not be like this in the future
CMD ["sleep","3600"]
