setenv('BLAS_VERSION','/usr/lib/libblas.so');
setenv('LAPACK_VERSION','/usr/lib/liblapack.so');
mex -largeArrayDims ../src/matlabinterface.cpp -L../libs/ -I../include/ -v -lfreewake -llapack ../src/matlabinterface.cpp
[surfs] = matlabinterface();
plotState(surfs)