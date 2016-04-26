#include <cmath>
#include <cstdlib>
#include "matlabinterface.h"

extern void _main();
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
        //Run Simulation
        SimulationManager sm = SimulationManager();
        sm.setReferenceVelocity( 100.0 );
        sm.setDt( 0.01 );
        sm.setReferenceSurface( ReferenceSurface( 10.0, 10.0, 1.0) );
        LiftingSurface ls = LiftingSurface(20,3,50);
        ls.setFreeWake( true );
        ls.setAspectRatio( 10.0 );
        ls.setPitch( 5.0 * M_PI / 180.0 );
        ls.getHorseshoeLattice().setHasTrailers(false); 
        ls.updateLattice( );
        sm.addSurface(&ls);
        sm.setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        int nStep = 0;
       // mexprintf("Solve\n");
        while ( nStep < 100 ){
            sm.step();
            nStep++; 
        }
        
        int nFields = 6;
        int nSurfaces = 1;
        const char* fieldNames[] = {"xSurface","ySurface","zSurface","xWake","yWake", "zWake"}; 
        plhs[0] = mxCreateStructMatrix(1, nSurfaces , nFields, fieldNames);
        for(int iSurface = 0; iSurface < nSurfaces; iSurface++){
            HorseshoeLattice& hl = ls.getHorseshoeLattice();  
            //Fill endpoint matrices
            int ni = ls.getHorseshoeLattice().ni()+1;
            int nj = ls.getHorseshoeLattice().nj()+1;
            int maxEndpoints = ni*nj; 
            mxArray* xmat = mxCreateDoubleMatrix(ni,nj, mxREAL);
            mxArray* ymat = mxCreateDoubleMatrix(ni,nj, mxREAL);
            mxArray* zmat = mxCreateDoubleMatrix(ni,nj, mxREAL);
            double* x = mxGetPr(xmat);
            double* y = mxGetPr(ymat);
            double* z = mxGetPr(zmat);
            for (int i = 0; i<ni; i++){ 
                for (int j = 0; j < nj; j++){
                    int n = j*ni+i;
                    x[n] = hl.getEndPoints()[i][j].x;
                    y[n] = hl.getEndPoints()[i][j].y;
                    z[n] = hl.getEndPoints()[i][j].z;
                }
            }
            mxSetFieldByNumber(plhs[0],iSurface,0,xmat);
            mxSetFieldByNumber(plhs[0],iSurface,1,ymat);
            mxSetFieldByNumber(plhs[0],iSurface,2,zmat);

            //Fill wake endpoint matrices
            VortexLattice& vl = ls.getVortexLattice();  
            ni = vl.ni();
            nj = vl.nj();
            maxEndpoints = ni*nj; 
            xmat = mxCreateDoubleMatrix(ni,nj, mxREAL);
            ymat = mxCreateDoubleMatrix(ni,nj, mxREAL);
            zmat = mxCreateDoubleMatrix(ni,nj, mxREAL);
            x = mxGetPr(xmat);
            y = mxGetPr(ymat);
            z = mxGetPr(zmat);
            for (int i = 0; i<ni; i++){ 
                for (int j = 0; j < nj; j++){
                    int n = j*ni+i;
                    x[n] = vl.endPoints()[i][j].x;
                    y[n] = vl.endPoints()[i][j].y;
                    z[n] = vl.endPoints()[i][j].z;
                }
            }
            mxSetFieldByNumber(plhs[0],iSurface,3,xmat);
            mxSetFieldByNumber(plhs[0],iSurface,4,ymat);
            mxSetFieldByNumber(plhs[0],iSurface,5,zmat);

        }
        //plhs[0] = xmat;
        //plhs[1] = ymat;
        //plhs[2] = zmat;
        
}

