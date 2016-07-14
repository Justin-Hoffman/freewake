#include <cmath>
#include <cstdlib>
#include <string.h>
#include "matlabinterface.h"

extern void _main();
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
        
        MatlabInterfaceStruct inArgs = validateArgs( nrhs, prhs);

        double omega = inArgs.omega;
        double dt = inArgs.dt;
        ReferenceSurface refS = ReferenceSurface( inArgs.refA, inArgs.refL, inArgs.refC);
        refS.pgCorrection = inArgs.doPrandtlGlauert;
        refS.vMach = inArgs.vMach;

        //Run Simulation
        SimulationManager sm = SimulationManager();
        sm.setReferenceVelocity( inArgs.refV );
        sm.setDt( dt );
        sm.setReferenceSurface( refS );
        LiftingSurface *lsp = new LiftingSurface(inArgs.nSpan,inArgs.nChord,inArgs.nNearWake,inArgs.nFarWake);
        std::vector<LiftingSurface> surfs = std::vector<LiftingSurface>(inArgs.nSurfaces, LiftingSurface(inArgs.nSpan,inArgs.nChord,inArgs.nNearWake,inArgs.nFarWake) );
        for (int i = 0; i < inArgs.nSurfaces; i++){
            LiftingSurface &ls = surfs[i];
            ls.setFreeWake( inArgs.isFreeWake );
            ls.setFreeTipVortex( inArgs.isFreeTipVortex );
            ls.setAspectRatio( inArgs.surfaceAR );
            ls.setPitch( inArgs.surfacePitch );
            ls.setCoreRadius( 1E-4 );
            ls.setTipDihedral( inArgs.surfaceTipDihedral );
            ls.setTipDihedralBreak( inArgs.surfaceTipDihedralBreak );
            ls.getHorseshoeLattice().spanwiseSpacing(PointSpacing::Cosine);
            ls.getHorseshoeLattice().chordwiseSpacing(PointSpacing::Cosine);
            ls.getHorseshoeLattice().setHasTrailers( inArgs.hasFixedTrailers ); 
            ls.updateLattice( );
            ls.getHorseshoeLattice().translate( Vec3D(0.0, 0.75, 0.0) );
            ls.getHorseshoeLattice().rotate(  Vec3D(0.0,0.0,0.0), Vec3D(0.0, 0.0, -1.0), ( (double)  i ) * 2.0 * M_PI/((double) inArgs.nSurfaces) );
            ls.getVortexLattice().fixToTrailingEdge( ls.getHorseshoeLattice() );
            ls.getVortexLattice().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.0 );
            ls.getVortexLattice().fixToTrailingEdge( ls.getHorseshoeLattice() );
            ls.getTipFilament().fixToWake( ls.getVortexLattice() );
            ls.getTipFilament().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.0 );
            sm.addSurface(&ls);
        }
        //LiftingSurface ls2 = LiftingSurface(ls);
        //ls2.setFreeWake( inArgs.isFreeWake );
        //ls2.setFreeTipVortex( inArgs.isFreeTipVortex );
        //ls2.getHorseshoeLattice().rotate(  Vec3D(0.0,0.0,0.0), Vec3D(0.0, 0.0, -1.0), M_PI );
        //ls2.getVortexLattice().fixToTrailingEdge( ls2.getHorseshoeLattice() );
        //ls2.getVortexLattice().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.0 );
        //ls2.getVortexLattice().fixToTrailingEdge( ls2.getHorseshoeLattice() );
        //ls2.getTipFilament().fixToWake( ls2.getVortexLattice() );
        //ls2.getTipFilament().initializeToHelix( Vec3D(0.0, 0.0, 1.0), (omega)*sm.dt(), -0.0 );
        //sm.addSurface(&ls2);

        sm.setGlobalLinearVelocity( inArgs.globalLinearVelocity.rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationAxis( inArgs.globalRotationAxis.rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -0.0*M_PI/180.0 ) );
        sm.setGlobalRotationRate( omega );
        
        int nt = inArgs.nt;

        int nFields = 24;
        int nSurfaces = sm.getNSurfaces();
        const char* fieldNames[] = {"xSurface","ySurface","zSurface","xWake","yWake", "zWake","xCp","yCp","zCp","xTipFilament","yTipFilament", "zTipFilament", "rcTip", "xSpanwiseForce", "ySpanwiseForce", "zSpanwiseForce", 
                                    "CFX","CFY","CFZ","CMX","CMY","CMZ", "T", "InputParameters"}; 
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
        mxArray** rctf  = (mxArray**) malloc( nSurfaces * sizeof( mxArray* ) );
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
            rctf[iSurface]= mxCreateDoubleMatrix(ni,nj-1, mxREAL);


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
            switch ( inArgs.integrationScheme ) {
                case IntegrationScheme::EULER : sm.step(); break;
                case IntegrationScheme::RK2   : sm.stepRK2(); break;
                case IntegrationScheme::PCC   : sm.stepPCC(); break;
                case IntegrationScheme::PC2B  : sm.stepPC2B(); break;
            }
            time += sm.dt();
                for(int iSurface = 0; iSurface < nSurfaces; iSurface++){
                    LiftingSurface &l = sm.getSurface( iSurface );
                    double* cfzp = mxGetPr( cfz[iSurface] );
                    double* cmzp = mxGetPr( cmz[iSurface] );
                    double* tp  = mxGetPr(  t[iSurface] );
                    cfzp[nStep]  = sm.forcesAndMoments().bodyForceCoeff.z;
                    cmzp[nStep]  = sm.forcesAndMoments().bodyMomentCoeff.z;
                    tp[nStep ]  = time;
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
                    double* rc= mxGetPr(rctf[iSurface] );
                    for (int i = 0; i<ni; i++){ 
                        for (int j = 0; j < nj; j++){
                            int n = j*ni+i;
                            x[n] = tf.endPoints()[i][j].x;
                            y[n] = tf.endPoints()[i][j].y;
                            z[n] = tf.endPoints()[i][j].z;
                            if ( j < nj-1 ){ 
                                rc[n] = tf.rc()[i][j];
                            }
                        }
                    }
                    mxSetFieldByNumber(plhs[0],iSurface,9 ,xtf[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,10,ytf[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,11,ztf[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,12,rctf[iSurface]);
                 
                    
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
                    mxSetFieldByNumber(plhs[0],iSurface,13,xspfc[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,14,yspfc[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,15,zspfc[iSurface]);

                    
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
                    mxSetFieldByNumber(plhs[0],iSurface,16, cfx[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,17, cfy[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,18, cfz[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,19, cmx[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,20, cmy[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,21, cmz[iSurface] );
                    mxSetFieldByNumber(plhs[0],iSurface,22, t[iSurface]);
                    mxSetFieldByNumber(plhs[0],iSurface,23, mxDuplicateArray(prhs[0]) );
                    
                } 
                mexCallMATLAB(0, NULL, 1, plhs, "plotState");
            }
            nStep++; 
        }
    free(xsurf); free(ysurf); free(zsurf);    
    free(xwake); free(ywake); free(zwake);    
    free(xcp); free(ycp); free(zcp);    
    free(xtf); free(ytf); free(ztf); free(rctf);
    free(xspfc); free(yspfc); free(zspfc);   
    free(cfx); free(cfy); free(cfz); 
    free(cmx); free(cmy); free(cmz); 
    free(t); 
}

MatlabInterfaceStruct validateArgs( int nrhs, const mxArray *prhs[] ) {
    MatlabInterfaceStruct inArgs;
    mxArray *val;
    if (nrhs != 1 ){
        mexErrMsgIdAndTxt("Freewake:InvalidArgumentCount", "Function called with other than 1 argument.");
    }
    
    val = mxGetField( prhs[0], 0, "omega");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: omega");
    } else {
        inArgs.omega = *mxGetPr(val);
    }
    
    val = mxGetField( prhs[0], 0, "dt");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: dt");
    } else {
        inArgs.dt = *mxGetPr(val);
    }
   
    val = mxGetField( prhs[0], 0, "nt");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: nt");
    } else {
        inArgs.nt = (int) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "globalLinearVelocity");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: globalLinearVelocity");
    } else {
        double* pr = mxGetPr(val);
        inArgs.globalLinearVelocity = Vec3D( pr[0], pr[1], pr[2]);
    }

    val = mxGetField( prhs[0], 0, "globalRotationAxis");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: globalRotationAxis");
    } else {
        double* pr = mxGetPr(val);
        inArgs.globalRotationAxis = Vec3D( pr[0], pr[1], pr[2] );
    }

    val = mxGetField( prhs[0], 0, "refL");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: refL");
    } else {
        inArgs.refL = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "refC");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: refC");
    } else {
        inArgs.refC = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "refA");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: refA");
    } else {
        inArgs.refA = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "refV");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: refV");
    } else {
        inArgs.refV = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "vMach");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: vMach");
    } else {
        inArgs.vMach = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "nSurfaces");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: nSurfaces");
    } else {
        inArgs.nSurfaces = (int) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "nChord");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: nChord");
    } else {
        inArgs.nChord = (int) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "nSpan");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: nSpan");
    } else {
        inArgs.nSpan = (int) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "nNearWake");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: nNearWake");
    } else {
        inArgs.nNearWake = (int) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "nFarWake");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: nFarWake");
    } else {
        inArgs.nFarWake = (int) *mxGetPr(val);
    } 

    val = mxGetField( prhs[0], 0, "surfaceAR");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: surfaceAR");
    } else {
        inArgs.surfaceAR = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "surfacePitch");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: surfacePitch");
    } else {
        inArgs.surfacePitch = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "surfaceTipDihedral");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: surfaceTipDihedral");
    } else {
        inArgs.surfaceTipDihedral = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "surfaceTipDihedralBreak");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: surfaceTipDihedralBreak");
    } else {
        inArgs.surfaceTipDihedralBreak = *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "isFreeWake");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: isFreeWake");
    } else {
        inArgs.isFreeWake = (bool) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "isFreeTipVortex");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: isFreeTipVortex");
    } else {
        inArgs.isFreeTipVortex = (bool) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "hasFixedTrailers");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: hasFixedTrailers");
    } else {
        inArgs.hasFixedTrailers = (bool) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "doPrandtlGlauert");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: doPrandtlGlauert");
    } else {
        inArgs.doPrandtlGlauert = (bool) *mxGetPr(val);
    }

    val = mxGetField( prhs[0], 0, "integrationScheme");
    if (val == 0){
        mexErrMsgIdAndTxt("Freewake:ErrorReadingField", "The following field is missing or incorrectly formatted: integrationScheme");
    } else {
        char* scheme = mxArrayToString(val);
        if ( !strcmp( scheme, "EULER") ){
            inArgs.integrationScheme = IntegrationScheme::EULER;
        } else if ( !strcmp( scheme, "RK2") ) {
            inArgs.integrationScheme = IntegrationScheme::RK2;
        } else if ( !strcmp( scheme, "PCC") ) {
            inArgs.integrationScheme = IntegrationScheme::PCC;
        } else if ( !strcmp( scheme, "PC2B") ) {
            inArgs.integrationScheme = IntegrationScheme::PC2B;
        }
    }


    return inArgs;
}

