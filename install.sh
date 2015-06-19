
cmake

-DCMAKE_C_COMPILER=cc
-DCMAKE_CXX_COMPILER=CC
-DCMAKE_Fortran_COMPILER=ftn
-DBLAS_LIB="/opt/cray/libsci/12.0.00/GNU/47/interlagos/lib/libsci_gnu_mp.a"
-DLAPACK_LIB=/opt/cray/libsci/12.0.00/GNU/47/interlagos/lib/libsci_gnu_mp.a
-DBUILD_SHARED_LIBS=OFF
-DCMAKE_EXE_LINKER_FLAGS="-static"
-DCMAKE_INSTALL_PREFIX=${HOME}/dealii_install
-DDEAL_II_WITH_PETSC=ON
-DMPI_LIBRARIES=${MPICH_DIR}/lib
-DMPI_LINKER_FLAGS=libmpichcxx_gnu_49_mt.a
-DMPI_INCLUDE_DIRS=${MPICH_DIR}/include   ..


mkdir build
cd build
cmake
-DCMAKE_C_COMPILER="gcc"
-DMPI_C_LIBRARIES="/opt/cray/mpt/default/gni/mpich2-GNU/48/lib"
-DMPI_C_INCLUDE_PATH="/opt/cray/mpt/default/gni/mpich2-GNU/48/include"
-DCMAKE_CXX_COMPILER="g++" \
-DMPI_CXX_LIBRARIES="/opt/cray/mpt/default/gni/mpich2-GNU/48/lib" \
-DMPI_CXX_INCLUDE_PATH="/opt/cray/mpt/defaultgni/mpich2-GNU/48/include" \
-DCMAKE_Fortran_COMPILER="gfortran" \
-DMPI_Fortran_LIBRARIES="/opt/cray/mpt/6.0.0/gni/mpich2-GNU/48/lib" \
-DMPI_Fortran_INCLUDE_PATH="/opt/cray/mpt/6.0.0/gni/mpich2-GNU/48/include" \
-DCMAKE_INSTALL_PREFIX=$HOME/dealII/install \
-DDEAL_II_WITH_THREADS=OFF \
-DDEAL_II_WITH_MPI=ON \
-DDEAL_II_WITH_UMFPACK=ON \
-DDEAL_II_WITH_LAPACK=ON \
-DLAPACK_DIR=/opt/intel/composerxe/mkl/lib/intel64 \
-DDEAL_II_WITH_TRILINOS=OFF \
-DDEAL_II_WITH_PETSC=ON \
-DDEAL_II_WITH_METIS=OFF \
-DDEAL_II_WITH_ZLIB=ON \
-DBUILD_SHARED_LIBS=ON \
-DDEAL_II_WITH_P4EST=ON \
../source