#include <cmath>
#include <cstdlib>
#include "matlabinterface.h"

extern void _main();
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
        double omega = 176.0; //Roughly mach 0.6 for a 7.5 ft span rotor (Caradonna Tung Rotor)
        double dt = 2.0*M_PI/ omega / 50.0;
        double r = 6.0; 
        double c = 1.0;
        double refA = M_PI * r * r;
        //Run Simulation
        SimulationManager sm = SimulationManager();
        sm.setReferenceVelocity( omega * r );
        sm.setDt( dt );
        sm.setReferenceSurface( ReferenceSurface( refA, r, c) );
        LiftingSurface ls = LiftingSurface(12,5,200);
        ls.setFreeWake( true );
        ls.setAspectRatio( 5.25 );
        ls.setPitch( 5.0 * M_PI / 180.0 );
        ls.getHorseshoeLattice().setHasTrailers(false); 
        ls.updateLattice( );
        ls.getHorseshoeLattice().translate( Vec3D(0.0, 0.75, 0.0) );
        ls.getVortexLattice().fixToTrailingEdge( ls.getHorseshoeLattice() );
        ls.getVortexLattice().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.03 );

        LiftingSurface ls2 = LiftingSurface(ls);
        ls2.setFreeWake( true );
        ls2.getHorseshoeLattice().rotate(  Vec3D(0.0,0.0,0.0), Vec3D(0.0, 0.0, -1.0), M_PI );
        ls2.getVortexLattice().fixToTrailingEdge( ls2.getHorseshoeLattice() );
        ls2.getVortexLattice().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.03 );
        
        
        sm.addSurface(&ls);
        sm.addSurface(&ls2);
        sm.setGlobalLinearVelocity( Vec3D(0.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationAxis( Vec3D(0.0, 0.0, 1.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationRate( omega );
        
        int nt = 1000;
        int nFields = 14;
        int nSurfaces = sm.getNSurfaces();
        const char* fieldNames[] = {"xSurface","ySurface","zSurface","xWake","yWake", "zWake","xCp","yCp","zCp", "xSpanwiseForce", "ySpanwiseForce", "zSpanwiseForce", "CT", "T"}; 
        plhs[0] = mxCreateStructMatrix(1, nSurfaces , nFields, fieldNames);
        mxArray** xsurf = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** ysurf = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** zsurf = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** xwake = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** ywake = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** zwake = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** xcp   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** ycp   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** zcp   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** xspfc = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** yspfc = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** zspfc = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        
        mxArray** ct    = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** t     = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        
        
        for(int iSurface = 0; iSurface < nSurfaces; iSurface++){
            LiftingSurface &l = sm.getSurface( iSurface );
            HorseshoeLattice &hl = l.getHorseshoeLattice();  
            //Fill endpoint matrices
            int ni = hl.ni()+1;
            int nj = hl.nj()+1;
            xsurf[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);
            ysurf[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);
            zsurf[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);

            VortexLattice& vl = l.getVortexLattice();  
            ni = vl.ni();
            nj = vl.nj();
            xwake[iSurface]  = mxCreateDoubleMatrix(ni,nj, mxREAL);
            ywake[iSurface]  = mxCreateDoubleMatrix(ni,nj, mxREAL);
            zwake[iSurface]  = mxCreateDoubleMatrix(ni,nj, mxREAL);
            
            ni = hl.ni();
            nj = hl.nj();
            xcp[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);
            ycp[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);
            zcp[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);

            ni = hl.ni();
            nj = hl.nj();
            xspfc[iSurface]  = mxCreateDoubleMatrix(ni, 1, mxREAL);
            yspfc[iSurface]  = mxCreateDoubleMatrix(ni, 1, mxREAL);
            zspfc[iSurface]  = mxCreateDoubleMatrix(ni, 1, mxREAL);
            

            ct[iSurface]     = mxCreateDoubleMatrix(nt, 1, mxREAL);
            t[iSurface]      = mxCreateDoubleMatrix(nt, 1, mxREAL);
        }

        int nStep = 0;
        double time = 0;
        // mexprintf("Solve\n");
        while ( nStep < nt ){
            sm.step();
            time += sm.dt();
            if (nStep % 1 == 0){
                for(int iSurface = 0; iSurface < nSurfaces; iSurface++){
                    LiftingSurface &l = sm.getSurface( iSurface );
                    HorseshoeLattice &hl = l.getHorseshoeLattice();  
                    //Fill endpoint matrices
                    int ni = hl.ni()+1;
                    int nj = hl.nj()+1;
                    int maxEndpoints = ni*nj; 
                    double* x = mxGetPr(xsurf[iSurface] );
                    double* y = mxGetPr(ysurf[iSurface] );
                    double* z = mxGetPr(zsurf[iSurface] );
                    for (int i = 0; i<ni; i++){ 
                        for (int j = 0; j < nj; j++){
                            int n = j*ni+i;
                            x[n] = hl.getEndPoints()[i][j].x;
                            y[n] = hl.getEndPoints()[i][j].y;
                            z[n] = hl.getEndPoints()[i][j].z;
                        }
                    }
                    mxSetFieldByNumber(plhs[0],iSurface,0,xsurf[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,1,ysurf[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,2,zsurf[iSurface]);

                    //Fill wake endpoint matrices
                    VortexLattice& vl = l.getVortexLattice();  
                    ni = vl.ni();
                    nj = vl.nj();
                    maxEndpoints = ni*nj; 
                    x = mxGetPr(xwake[iSurface] );
                    y = mxGetPr(ywake[iSurface] );
                    z = mxGetPr(zwake[iSurface] );
                    for (int i = 0; i<ni; i++){ 
                        for (int j = 0; j < nj; j++){
                            int n = j*ni+i;
                            x[n] = vl.endPoints()[i][j].x;
                            y[n] = vl.endPoints()[i][j].y;
                            z[n] = vl.endPoints()[i][j].z;
                        }
                    }
                    mxSetFieldByNumber(plhs[0],iSurface,3,xwake[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,4,ywake[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,5,zwake[iSurface]);
                    
                    //Fill control point matrices
                    ni = hl.ni();
                    nj = hl.nj();
                    maxEndpoints = ni*nj; 
                    x = mxGetPr(xcp[iSurface] );
                    y = mxGetPr(ycp[iSurface] );
                    z = mxGetPr(zcp[iSurface] );
                    for (int i = 0; i<ni; i++){ 
                        for (int j = 0; j < nj; j++){
                            int n = j*ni+i;
                            x[n] = hl.getControlPoints()[i][j].x;
                            y[n] = hl.getControlPoints()[i][j].y;
                            z[n] = hl.getControlPoints()[i][j].z;
                        }
                    }
                    mxSetFieldByNumber(plhs[0],iSurface,6,xcp[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,7,ycp[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,8,zcp[iSurface]);

                    
                    //Fill spanwise force matrices
                    ni = vl.ni();
                    nj = vl.nj();
                    x = mxGetPr(xspfc[iSurface] );
                    y = mxGetPr(yspfc[iSurface] );
                    z = mxGetPr(zspfc[iSurface] );
                    for (int i = 0; i<ni; i++){ 
                        x[i] = l.spanwiseForce()[i].x;
                        y[i] = l.spanwiseForce()[i].y;
                        z[i] = l.spanwiseForce()[i].z;
                    }
                    mxSetFieldByNumber(plhs[0],iSurface,9,xspfc[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,10,yspfc[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,11,zspfc[iSurface]);

                    
                    double* ctp = mxGetPr( ct[iSurface] );
                    double* tp  = mxGetPr(  t[iSurface] );
                    ctp[nStep]  = sm.forcesAndMoments().bodyForceCoeff.z;
                    tp[nStep ]  = time;
                    mxSetFieldByNumber(plhs[0],iSurface,12,ct[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,13, t[iSurface]);
                    
                } 
                mexCallMATLAB(0, NULL, 1, plhs, "plotState");
            }
            nStep++; 
        }
    free(xsurf); free(ysurf); free(zsurf);    
    free(xwake); free(ywake); free(zwake);    
    free(xcp); free(ycp); free(zcp);    
    free(xspfc); free(yspfc); free(zspfc);   
    free(ct); free(t); 
}

