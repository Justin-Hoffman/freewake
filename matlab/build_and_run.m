setenv('BLAS_VERSION','/usr/lib/libblas.so');
setenv('LAPACK_VERSION','/usr/lib/liblapack.so');
mex -largeArrayDims -L../libs/ -I../include/ -v -lfreewake -llapack -llapacke ../src/matlabinterface.cpp
[surfs] = matlabinterface();
plotState(surfs)