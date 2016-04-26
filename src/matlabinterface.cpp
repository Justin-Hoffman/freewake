#include <cmath>
#include <cstdlib>
#include "matlabinterface.h"

extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
        //Run Simulation
        SimulationManager sm = SimulationManager();
        sm.setReferenceVelocity( 100.0 );
        sm.setDt( 0.1 );
        sm.setReferenceSurface( ReferenceSurface( 10.0, 10.0, 1.0) );
        LiftingSurface ls = LiftingSurface(10,1,10);
        ls.setFreeWake( true );
        ls.setAspectRatio( 10.0 );
        ls.setPitch( 0.0 * M_PI / 180.0 );
        ls.getHorseshoeLattice().setHasTrailers(false); 
        HorseshoeLattice& hl = ls.getHorseshoeLattice();  
        ls.updateLattice( );
        sm.addSurface(&ls);
        sm.setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -10.0*M_PI/180.0 ) );
        int nStep = 0;
        /*printf("Solve\n");
        while ( nStep < 100 ){
            sm.step();
            nStep++; 
        }
        
        int nFields = 3; 
        mxClassID* classIDFlags = (mxClassID*) mxCalloc( nFields, sizeof(mxClassID) );
        for( int iField = 0; iField < nFields; iField++){
            classIDFlags[iField] = mxDOUBLE_CLASS; 
        }
       
        int ni = ls.getHorseshoeLattice().ni()+1;
        int nj = ls.getHorseshoeLattice().nj()+1;

        int maxEndpoints = ni*nj; 

        mxArray* xmat = mxCreateNumericMatrix(ni,nj, mxDOUBLE_CLASS, mxREAL);
        mxArray* ymat = mxCreateNumericMatrix(ni,nj, mxDOUBLE_CLASS, mxREAL);
        mxArray* zmat = mxCreateNumericMatrix(ni,nj, mxDOUBLE_CLASS, mxREAL);
        for (int i = 0; i<ni; i++){ 
            for (int j = 0; j < nj; j++){
                double *x = mxGetPr(xmat);
                double *y = mxGetPr(ymat);
                double *z = mxGetPr(zmat);
                int n = i*ni+j;
                x[n] = hl.getEndPoints()[i][j].x;
                y[n] = hl.getEndPoints()[i][j].y;
                z[n] = hl.getEndPoints()[i][j].z;
            }
        }
        plhs[0] = xmat;
        plhs[1] = ymat;
        plhs[2] = zmat;
        */
}

