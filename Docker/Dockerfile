# Base HPC SDK image
FROM nvcr.io/nvidia/nvhpc:24.3-devel-cuda_multi-ubuntu22.04

# Set noninteractive frontend
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build dependencies for OpenFOAM
RUN apt-get update && apt-get install -y --no-install-recommends \
    git ca-certificates wget software-properties-common \
    build-essential cmake python3 python3-dev \
    flex bison libfl-dev \
    zlib1g-dev libreadline-dev libncurses-dev \
    libopenmpi-dev openmpi-bin \
    libscotch-dev libscotchmetis-dev \
    && rm -rf /var/lib/apt/lists/*


# Enable HPC SDK environment
ENV PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/bin:$PATH
ENV PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/comm_libs/mpi/bin:$PATH
ENV MANPATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/man:$MANPATH
ENV LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/lib:$LD_LIBRARY_PATH
ENV CC=nvc
ENV CXX=nvc++
ENV FC=nvfortran

# Clone OpenFOAM (HTTPS fallback; replace with SSH if you use your keys)
WORKDIR /opt/OpenFOAM
# RUN git clone --branch of-stdpar https://github.com/alpha-unito/OpenFOAM.git
COPY . .

# Set OpenFOAM environment
ENV WM_PROJECT_DIR=/opt/OpenFOAM
ENV FOAM_INST_DIR=/opt
ENV WM_COMPILER=Clang
ENV WM_MPLIB=SYSTEMMPI
ENV GPU_PLATFORMS=cc90,cc86

# Source environment for building
RUN bash -c "source /opt/OpenFOAM/etc/bashrc && \
    foamSystemCheck"

# Build OpenFOAM
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc && \
    cd $WM_PROJECT_DIR/wmake/src && \
    make "


RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc && \
    # OS-specific
    cd $WM_PROJECT_DIR/src/OSspecific/POSIX && ./Allwmake -j 4 && \
    cd $WM_PROJECT_DIR/src/Pstream && ./Allwmake -j 8 "


    # Core OpenFOAM libraries
# OpenFOAM core
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/OpenFOAM && wmake -j 8 ."

# fileFormats
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/fileFormats && wmake -j 8 ."

# surfMesh
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/surfMesh && wmake -j 8 ."

# meshTools
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/meshTools && wmake -j 8 ."

# finiteArea
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/finiteArea && wmake -j 8 ."

# finiteVolume
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/finiteVolume && wmake -j 8 ."

# fused/finiteVolume
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/fused/finiteVolume && wmake -j 8 ."

# mesh/blockMesh
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/mesh/blockMesh && wmake -j 8 ."

# mesh/extrudeModel  (requires blockMesh)
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/mesh/extrudeModel && wmake -j 8 ."

# dynamicMesh  (requires extrudeModel)
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/dynamicMesh && wmake -j 8 ."


# dynamicMesh  (requires setFields)
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc \
    && cd $WM_PROJECT_DIR/src/genericPatchFields && wmake -j 8 ."    

# Build OpenFOAM application blockMesh
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc && \
    cd $WM_PROJECT_DIR/applications/utilities/mesh/generation/blockMesh && \
    wmake "

# Build OpenFOAM application setField
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc && \
    cd $WM_PROJECT_DIR/applications/utilities/preProcessing/setFields && \
    wmake "

# Build OpenFOAM application laplacianFoam
RUN bash -c "source $WM_PROJECT_DIR/etc/bashrc && \
    cd $WM_PROJECT_DIR/TestLaplacianFoam/laplacianFoam && \
    wmake "


# Set environment on container startup
RUN echo 'source /opt/OpenFOAM/etc/bashrc' >> /etc/bash.bashrc

WORKDIR /root
CMD ["/bin/bash"]
