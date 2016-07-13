setenv('BLAS_VERSION','/usr/lib/libblas.so');
setenv('LAPACK_VERSION','/usr/lib/liblapack.so');
%mex -largeArrayDims ../src/matlabinterface.cpp -L../libs/ -I../include/ -v -lfreewake -llapack -llapacke ../src/matlabinterface.cpp
mex -largeArrayDims -L../libs/ -I../include/ -v -lfreewake -llapack -llapacke ../src/matlabinterface.cpp

inArgs = defaultRun();
nWake = 100:50:250;
for i = 1:length(nWake)
    inArgs.nNearWake = nWake(i);
    [surfs] = matlabinterface(inArgs);
    fname = ['PCC-' num2str(nWake(i)) '-50-per-rev.mat'];
    save(fname,'surfs');
end
plotState(surfs)