setenv('BLAS_VERSION','/usr/lib/libblas.so');
setenv('LAPACK_VERSION','/usr/lib/liblapack.so');
mex -largeArrayDims ../src/matlabinterface.cpp -L../libs/ -I../include/ -v -lfreewake -llapack -llapacke ../src/matlabinterface.cpp
%mex -largeArrayDims -L../libs/ -I../include/ -v -lfreewake -llapack -llapacke ../src/matlabinterface.cpp

inArgs = defaultRun();
nWake = 50:50:50;
tipCant = 0:15:0;
for i = 1:length(nWake)
    for j = 1:length(tipCant)
        inArgs.nNearWake = nWake(i);
        [surfs] = matlabinterface(inArgs);
        fname = [inArgs.integrationScheme '-' num2str(nWake(i)) '-' num2str(tipCant(j)) 'cant-50-per-rev(B).mat'];
        save(fname,'surfs');
    end
end
plotState(surfs)