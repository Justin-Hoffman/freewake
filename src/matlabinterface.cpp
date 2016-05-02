#include <cmath>
#include <cstdlib>
#include "matlabinterface.h"

extern void _main();
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
        double omega = 176.0; //Roughly mach 0.6 for a 7.5 ft span rotor (Caradonna Tung Rotor)
        double dt = 2.0*M_PI/ omega / 72.0;
        double r = 6.0; 
        double c = 1.0;
        double refA = M_PI * r * r;
        //Run Simulation
        SimulationManager sm = SimulationManager();
        sm.setReferenceVelocity( omega * r );
        sm.setDt( dt );
        sm.setReferenceSurface( ReferenceSurface( refA, r, c) );
        LiftingSurface ls = LiftingSurface(20,10,6,250);
        ls.setFreeWake( true );
        ls.setFreeTipVortex( true );
        ls.setAspectRatio( 5.25 );
        ls.setPitch( 5.0 * M_PI / 180.0 );
        ls.setCoreRadius( 1E-3 );
        ls.setTipDihedral( 0.0 * M_PI/180.0 );
        ls.setTipDihedralBreak( .9*5.25/6.0 );
        ls.getHorseshoeLattice().spanwiseSpacing(PointSpacing::Cosine);
        ls.getHorseshoeLattice().chordwiseSpacing(PointSpacing::Cosine);
        ls.getHorseshoeLattice().setHasTrailers(false); 
        ls.updateLattice( );
        ls.getHorseshoeLattice().translate( Vec3D(0.0, 0.75, 0.0) );
        ls.getVortexLattice().fixToTrailingEdge( ls.getHorseshoeLattice() );
        ls.getVortexLattice().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.0 );
        ls.getVortexLattice().fixToTrailingEdge( ls.getHorseshoeLattice() );
        ls.getTipFilament().fixToWake( ls.getVortexLattice() );
        ls.getTipFilament().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.01 );

        LiftingSurface ls2 = LiftingSurface(ls);
        ls2.setFreeWake( true );
        ls.setFreeTipVortex( true );
        ls2.getHorseshoeLattice().rotate(  Vec3D(0.0,0.0,0.0), Vec3D(0.0, 0.0, -1.0), M_PI );
        ls2.getVortexLattice().fixToTrailingEdge( ls2.getHorseshoeLattice() );
        ls2.getVortexLattice().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.0 );
        ls2.getVortexLattice().fixToTrailingEdge( ls2.getHorseshoeLattice() );
        ls2.getTipFilament().fixToWake( ls2.getVortexLattice() );
        ls2.getTipFilament().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.01);
        
        
        sm.addSurface(&ls);
        sm.addSurface(&ls2);
        sm.setGlobalLinearVelocity( Vec3D(0.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationAxis( Vec3D(0.0, 0.0, 1.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationRate( omega );
        
        int nt = 2000;
        int nFields = 22;
        int nSurfaces = sm.getNSurfaces();
        const char* fieldNames[] = {"xSurface","ySurface","zSurface","xWake","yWake", "zWake","xCp","yCp","zCp","xTipFilament","yTipFilament", "zTipFilament", "xSpanwiseForce", "ySpanwiseForce", "zSpanwiseForce", 
                                    "CFX","CFY","CFZ","CMX","CMY","CMZ", "T"}; 
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
        mxArray** xtf   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** ytf   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** ztf   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** xspfc = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** yspfc = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** zspfc = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        
        mxArray** cfx   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** cfy   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** cfz   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** cmx   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** cmy   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** cmz   = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        mxArray** t     = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
        
        
        for(int iSurface = 0; iSurface < nSurfaces; iSurface++){
            LiftingSurface &l = sm.getSurface( iSurface );
            HorseshoeLattice &hl = l.getHorseshoeLattice();  
            TipFilament &tf      = l.getTipFilament();  
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
            
            ni = tf.ni();
            nj = tf.nj();
            xtf[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);
            ytf[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);
            ztf[iSurface] = mxCreateDoubleMatrix(ni,nj, mxREAL);


            ni = hl.ni();
            nj = hl.nj();
            xspfc[iSurface]  = mxCreateDoubleMatrix(ni, 1, mxREAL);
            yspfc[iSurface]  = mxCreateDoubleMatrix(ni, 1, mxREAL);
            zspfc[iSurface]  = mxCreateDoubleMatrix(ni, 1, mxREAL);
            

            cfx[iSurface]    = mxCreateDoubleMatrix(nt, 1, mxREAL);
            cfy[iSurface]    = mxCreateDoubleMatrix(nt, 1, mxREAL);
            cfz[iSurface]    = mxCreateDoubleMatrix(nt, 1, mxREAL);
            cmx[iSurface]    = mxCreateDoubleMatrix(nt, 1, mxREAL);
            cmy[iSurface]    = mxCreateDoubleMatrix(nt, 1, mxREAL);
            cmz[iSurface]    = mxCreateDoubleMatrix(nt, 1, mxREAL);
            t[iSurface]      = mxCreateDoubleMatrix(nt, 1, mxREAL);
        }

        int nStep = 0;
        double time = 0;
        // mexprintf("Solve\n");
        while ( nStep < nt ){
            sm.stepPC2B();
            time += sm.dt();
            if (nStep % 1 == 0){
                for(int iSurface = 0; iSurface < nSurfaces; iSurface++){
                    LiftingSurface &l = sm.getSurface( iSurface );
                    double* cfzp = mxGetPr( cfz[iSurface] );
                    double* cmzp = mxGetPr( cmz[iSurface] );
                    double* tp  = mxGetPr(  t[iSurface] );
                    cfzp[nStep]  = sm.forcesAndMoments().bodyForceCoeff.z;
                    cmzp[nStep]  = sm.forcesAndMoments().bodyMomentCoeff.z;
                    tp[nStep ]  = time;
                }
            }
            if (nStep % 5 == 0){
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
                   
                    //Fill tip vortex filament endpoints
                    TipFilament &tf      = l.getTipFilament();  
                    ni = tf.ni();
                    nj = tf.nj();
                    x = mxGetPr(xtf[iSurface] );
                    y = mxGetPr(ytf[iSurface] );
                    z = mxGetPr(ztf[iSurface] );
                    for (int i = 0; i<ni; i++){ 
                        for (int j = 0; j < nj; j++){
                            int n = j*ni+i;
                            x[n] = tf.endPoints()[i][j].x;
                            y[n] = tf.endPoints()[i][j].y;
                            z[n] = tf.endPoints()[i][j].z;
                        }
                    }
                    mxSetFieldByNumber(plhs[0],iSurface,9 ,xtf[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,10,ytf[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,11,ztf[iSurface]);
                 
                    
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
                    mxSetFieldByNumber(plhs[0],iSurface,12,xspfc[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,13,yspfc[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,14,zspfc[iSurface]);

                    
                    double* cfxp = mxGetPr( cfx[iSurface] );
                    double* cfyp = mxGetPr( cfy[iSurface] );
                    double* cfzp = mxGetPr( cfz[iSurface] );
                    double* cmxp = mxGetPr( cmx[iSurface] );
                    double* cmyp = mxGetPr( cmy[iSurface] );
                    double* cmzp = mxGetPr( cmz[iSurface] );
                    double* tp  = mxGetPr(  t[iSurface] );
                    cfxp[nStep]  = sm.forcesAndMoments().bodyForceCoeff.x;
                    cfyp[nStep]  = sm.forcesAndMoments().bodyForceCoeff.y;
                    cfzp[nStep]  = sm.forcesAndMoments().bodyForceCoeff.z;
                    cmxp[nStep]  = sm.forcesAndMoments().bodyMomentCoeff.x;
                    cmyp[nStep]  = sm.forcesAndMoments().bodyMomentCoeff.y;
                    cmzp[nStep]  = sm.forcesAndMoments().bodyMomentCoeff.z;
                    tp[nStep ]  = time;
                    mxSetFieldByNumber(plhs[0],iSurface,15, cfx[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,16, cfy[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,17, cfz[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,18, cmz[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,19, cmz[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,20, cmz[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,21, t[iSurface]);
                    
                } 
                mexCallMATLAB(0, NULL, 1, plhs, "plotState");
            }
            nStep++; 
        }
    free(xsurf); free(ysurf); free(zsurf);    
    free(xwake); free(ywake); free(zwake);    
    free(xcp); free(ycp); free(zcp);    
    free(xtf); free(ytf); free(ztf);    
    free(xspfc); free(yspfc); free(zspfc);   
    free(cfx); free(cfy); free(cfz); 
    free(cmx); free(cmy); free(cmz); 
    free(t); 
}

