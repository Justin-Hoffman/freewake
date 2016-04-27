#include <cmath>
#include <cstdlib>
#include "matlabinterface.h"

extern void _main();
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
        //Run Simulation
        SimulationManager sm = SimulationManager();
        sm.setReferenceVelocity( 1.0 );
        sm.setDt( 0.02 );
        sm.setReferenceSurface( ReferenceSurface( 10.0, 10.0, 1.0) );
        LiftingSurface ls = LiftingSurface(6,4,80);
        ls.setFreeWake( true );
        ls.setAspectRatio( 4.75 );
        ls.setPitch( 5.0 * M_PI / 180.0 );
        ls.getHorseshoeLattice().setHasTrailers(false); 
        ls.updateLattice( );
        ls.getHorseshoeLattice().translate( Vec3D(0.0, 1.25, 0.0) );
        LiftingSurface ls2 = LiftingSurface(ls);
        ls2.setFreeWake( true );
        ls2.getHorseshoeLattice().rotate(  Vec3D(0.0,0.0,0.0), Vec3D(0.0, 0.0, -1.0), M_PI );
        
        sm.addSurface(&ls);
        sm.addSurface(&ls2);
        sm.setGlobalLinearVelocity( Vec3D(0.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationAxis( Vec3D(0.0, 0.0, 1.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationRate( 2.0*M_PI );
        int nStep = 0;
       // mexprintf("Solve\n");
        while ( nStep < 250 ){
            sm.step();
            nStep++; 
        }
        
        int nFields = 6;
        int nSurfaces = sm.getNSurfaces();
        const char* fieldNames[] = {"xSurface","ySurface","zSurface","xWake","yWake", "zWake"}; 
        plhs[0] = mxCreateStructMatrix(1, nSurfaces , nFields, fieldNames);
        for(int iSurface = 0; iSurface < nSurfaces; iSurface++){
            LiftingSurface &l = sm.getSurface( iSurface );
            HorseshoeLattice &hl = l.getHorseshoeLattice();  
            //Fill endpoint matrices
            int ni = hl.ni()+1;
            int nj = hl.nj()+1;
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
            VortexLattice& vl = l.getVortexLattice();  
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

