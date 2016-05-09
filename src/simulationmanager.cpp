#include <cmath>
#include "simulationmanager.h"
#include "string.h"
#include "lapacke.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);

SimulationManager::SimulationManager() : needsSolve_( true ), surfaces_(), oldSurfaces_(), olderSurfaces_(), nOffset_(1,0), lastGamma_(), thisGamma_(), 
                                         refSurf_( ReferenceSurface(1.0, 1.0, 1.0) ), refV_(1.0), dt_( 0.10 ), 
                                         globalLinearVelocity_(), globalRotationAxis_(0.0,0.0,1.0), globalRotationRate_(0.0), 
                                         fomo_() {

}

SimulationManager::SimulationManager( const SimulationManager& s ) : needsSolve_( s.needsSolve_ ), surfaces_( s.surfaces_ ), oldSurfaces_( s.oldSurfaces_ ), olderSurfaces_( s.olderSurfaces_ ),  
                                                               nOffset_(s.nOffset_), lastGamma_( s.lastGamma_ ), thisGamma_( s.thisGamma_ ),
                                                               refSurf_( s.refSurf_ ), refV_( s.refV_ ), dt_( s.dt_ ), 
                                                               globalLinearVelocity_( s.globalLinearVelocity_), globalRotationAxis_( s.globalRotationAxis_ ), globalRotationRate_( s.globalRotationRate_),
                                                               fomo_( s.fomo_ ) {

}

SimulationManager::~SimulationManager(){
    for (int h = 0; h < oldSurfaces_.size(); h++){
        delete( oldSurfaces_[h] );
        delete( olderSurfaces_[h] );
    }
}

extern void printMatrix( char const* desc, int m, int n, double* a, int lda );
extern void print_int_vector( char const* desc, int n, int* a );

void SimulationManager::addSurface( LiftingSurface* s){
    surfaces_.push_back( s );
    oldSurfaces_.push_back( new LiftingSurface( *s ) );
    olderSurfaces_.push_back( new LiftingSurface( *s ) );
    nOffset_.push_back( s->getHorseshoeLattice().maxN() + nOffset_.back() );
    lastGamma_.push_back( std::vector<double>(s->nSpan()) ); //Make a double array of nSpan
    thisGamma_.push_back( std::vector<double>(s->nSpan()) ); //Make a double array of nSpan
}

void SimulationManager::setGlobalLinearVelocity( Vec3D v ){
    globalLinearVelocity_ = v;
}

void SimulationManager::setGlobalRotationAxis( Vec3D v ){
    globalRotationAxis_ = v;
}

void SimulationManager::setGlobalRotationRate( double r ){
    globalRotationRate_ = r;
}

Vec3D SimulationManager::getGlobalRotationAxis(){
    return globalRotationAxis_;
}

double SimulationManager::getGlobalRotationRate( ){
   return globalRotationRate_;
}

Vec3D SimulationManager::getGlobalLinearVelocity(){
    return globalLinearVelocity_;
}

void SimulationManager::step(){
    double dt = dt_;
    std::vector< LiftingSurface > originalSurfaces = std::vector< LiftingSurface >();
    for (int h = 0; h < (int) surfaces_.size(); h++){
        originalSurfaces.emplace_back( *surfaces_[h] );
    }
    dt_ = dt/2.0;

    //Do half time step for first part of RK2::
    solve();
    calculateWakeVelocities();
    advectWake();
    matchFirstWakePointsToSurface();
    calculateWakeVelocities();
    //Reset X,Y,Z,Gamma, but use new wake velocity
    for (int h = 0; h < (int) surfaces_.size(); h++){
        VortexLattice &vl = surfaces_[h]->getVortexLattice();
        VortexLattice &vlold = originalSurfaces[h].getVortexLattice();
        TipFilament &tv = surfaces_[h]->getTipFilament();
        TipFilament &tvold = originalSurfaces[h].getTipFilament();
        if (surfaces_[h]->freeWake()){
            vl.endPoints() = std::vector<std::vector<Vec3D> >(vlold.endPoints());
            vl.gammaI() = std::vector<std::vector< double > >(vlold.gammaI());
            vl.gammaJ() = std::vector<std::vector< double > >(vlold.gammaJ());
            vl.rcI() = std::vector<std::vector< double > >(vlold.rcI());
            vl.rcJ() = std::vector<std::vector< double > >(vlold.rcJ());
            tv.endPoints() = std::vector<std::vector<Vec3D> >(tvold.endPoints());
            tv.gamma() = std::vector<std::vector<double> >(tvold.gamma());
            tv.rc() = std::vector<std::vector<double> >(tvold.rc());
        }
    }
    double eps = 0.95;
    dt_ = dt;
    advectWake();
    fillWakeBC();
    
    //Apply explicit relaxation
    for (int h = 0; h < (int) surfaces_.size(); h++){
        VortexLattice &vl = surfaces_[h]->getVortexLattice();
        VortexLattice &vlold = originalSurfaces[h].getVortexLattice();
        TipFilament &tv = surfaces_[h]->getTipFilament();
        TipFilament &tvold = originalSurfaces[h].getTipFilament();
        if (surfaces_[h]->freeWake()){
            for( int i = 0; i < vl.ni(); i++ ){
                for( int j = 0; j < vl.nj(); j++ ){
                    vl.endPoints()[i][j] = (eps) * vl.endPoints()[i][j] + (1.0-eps) * vlold.endPoints()[i][j];
                    if ( j  < (vl.nj()-1) ) vl.gammaJ()[i][j] = (eps) * vl.gammaJ()[i][j] + (1.0-eps) * vlold.gammaJ()[i][j];
                }
            }
            for(int i = 0; i < 2; i++){
                for(int j = 0; j< tv.nj(); j++){
                    tv.endPoints()[i][j] = (eps) * tv.endPoints()[i][j] + (1.0-eps) * tvold.endPoints()[i][j];
                    if ( j  < (tv.nj()-1) ) tv.gamma()[i][j] = (eps) * tv.gamma()[i][j] + (1.0-eps) * tvold.gamma()[i][j];
                }
            }
        }
    }
    matchFirstWakePointsToSurface();
}

void SimulationManager::stepPCC(){
    std::vector< LiftingSurface > originalSurfaces = std::vector< LiftingSurface >();
    for (int h = 0; h < (int) surfaces_.size(); h++){
        originalSurfaces.emplace_back( *surfaces_[h] );
    }

    //Do time step to get velocities for PCC
    solve();
    calculateWakeVelocities();
    advectWake();
    matchFirstWakePointsToSurface();
    calculateWakeVelocities();
    //Reset X,Y,Z,Gamma, but use average wake velocity
    for (int h = 0; h < (int) surfaces_.size(); h++){
        VortexLattice &vl = surfaces_[h]->getVortexLattice();
        VortexLattice &vlold = originalSurfaces[h].getVortexLattice();
        TipFilament &tv = surfaces_[h]->getTipFilament();
        TipFilament &tvold = originalSurfaces[h].getTipFilament();
        if (surfaces_[h]->freeWake()){
            vl.endPoints() = std::vector<std::vector<Vec3D> >(vlold.endPoints());
            vl.gammaI() = std::vector<std::vector< double > >(vlold.gammaI());
            vl.gammaJ() = std::vector<std::vector< double > >(vlold.gammaJ());
            vl.rcI() = std::vector<std::vector< double > >(vlold.rcI());
            vl.rcJ() = std::vector<std::vector< double > >(vlold.rcJ());
            tv.endPoints() = std::vector<std::vector<Vec3D> >(tvold.endPoints());
            tv.gamma() = std::vector<std::vector<double> >(tvold.gamma());
            tv.rc() = std::vector<std::vector<double> >(tvold.rc());
            for(int i = 0; i < vl.ni(); i++){
                for( int j = 0; j < vl.nj(); j++){
                    vl.endPointVelocity()[i][j] = (vl.endPointVelocity()[i][j] + vlold.endPointVelocity()[i][j] ) /2.0;
                }
            }
            for(int i = 0; i < tv.ni(); i++){
                for( int j = 0; j < tv.nj(); j++){
                    tv.endPointVelocity()[i][j] = (tv.endPointVelocity()[i][j] + tvold.endPointVelocity()[i][j] ) /2.0;
                }
            }
        }
    }
    //Apply explicit relaxation
    advectWake();
    fillWakeBC();
    
    //Apply explicit relaxation
    double eps = 0.99;
    for (int h = 0; h < (int) surfaces_.size(); h++){
        VortexLattice &vl = surfaces_[h]->getVortexLattice();
        VortexLattice &vlold = originalSurfaces[h].getVortexLattice();
        TipFilament &tv = surfaces_[h]->getTipFilament();
        TipFilament &tvold = originalSurfaces[h].getTipFilament();
        if (surfaces_[h]->freeWake()){
            for( int i = 0; i < vl.ni(); i++ ){
                for( int j = 0; j < vl.nj(); j++ ){
                    vl.endPoints()[i][j] = (eps) * vl.endPoints()[i][j] + (1.0-eps) * vlold.endPoints()[i][j];
                    if ( j  < (vl.nj()-1) ) vl.gammaJ()[i][j] = (eps) * vl.gammaJ()[i][j] + (1.0-eps) * vlold.gammaJ()[i][j];
                }
            }
            for(int i = 0; i < 2; i++){
                for(int j = 0; j< tv.nj(); j++){
                    tv.endPoints()[i][j] = (eps) * tv.endPoints()[i][j] + (1.0-eps) * tvold.endPoints()[i][j];
                    if ( j  < (tv.nj()-1) ) tv.gamma()[i][j] = (eps) * tv.gamma()[i][j] + (1.0-eps) * tvold.gamma()[i][j];
                }
            }
        }
    }
    matchFirstWakePointsToSurface();
}
 
void SimulationManager::stepPC2B(){
    std::vector< LiftingSurface > originalSurfaces = std::vector< LiftingSurface >();
    for (int h = 0; h < (int) surfaces_.size(); h++){
        originalSurfaces.emplace_back( *surfaces_[h] );
    }
    //Do time step to get velocities for PCC
    solve();
    calculateWakeVelocities();
    advectWakePC2B();
    matchFirstWakePointsToSurface();
    calculateWakeVelocities();
    //Reset X,Y,Z,Gamma, but use average wake velocity
    for (int h = 0; h < (int) surfaces_.size(); h++){
        VortexLattice &vl = surfaces_[h]->getVortexLattice();
        VortexLattice &vlold = originalSurfaces[h].getVortexLattice();
        TipFilament &tv = surfaces_[h]->getTipFilament();
        TipFilament &tvold = originalSurfaces[h].getTipFilament();
        if (surfaces_[h]->freeWake()){
            vl.endPoints() = std::vector<std::vector<Vec3D> >(vlold.endPoints());
            vl.gammaI() = std::vector<std::vector< double > >(vlold.gammaI());
            vl.gammaJ() = std::vector<std::vector< double > >(vlold.gammaJ());
            vl.rcI() = std::vector<std::vector< double > >(vlold.rcI());
            vl.rcJ() = std::vector<std::vector< double > >(vlold.rcJ());
            tv.endPoints() = std::vector<std::vector<Vec3D> >(tvold.endPoints());
            tv.gamma() = std::vector<std::vector<double> >(tvold.gamma());
            tv.rc() = std::vector<std::vector<double> >(tvold.rc());
            for(int i = 0; i < vl.ni(); i++){
                for( int j = 0; j < vl.nj(); j++){
                    vl.endPointVelocity()[i][j] = (vl.endPointVelocity()[i][j] + vlold.endPointVelocity()[i][j] ) /2.0;
                }
            }
            for(int i = 0; i < tv.ni(); i++){
                for( int j = 0; j < tv.nj(); j++){
                    tv.endPointVelocity()[i][j] = (tv.endPointVelocity()[i][j] + tvold.endPointVelocity()[i][j] ) /2.0;
                }
            }
        }
    }
    advectWakePC2B();
    fillWakeBC();
    //Apply explicit relaxation
    double eps = 0.99;
    for (int h = 0; h < (int) surfaces_.size(); h++){
        VortexLattice &vl = surfaces_[h]->getVortexLattice();
        VortexLattice &vlold = originalSurfaces[h].getVortexLattice();
        TipFilament &tv = surfaces_[h]->getTipFilament();
        TipFilament &tvold = originalSurfaces[h].getTipFilament();
        if (surfaces_[h]->freeWake()){
            for( int i = 0; i < vl.ni(); i++ ){
                for( int j = 0; j < vl.nj(); j++ ){
                    vl.endPoints()[i][j] = (eps) * vl.endPoints()[i][j] + (1.0-eps) * vlold.endPoints()[i][j];
                    if ( j  < (vl.nj()-1) ) vl.gammaJ()[i][j] = (eps) * vl.gammaJ()[i][j] + (1.0-eps) * vlold.gammaJ()[i][j];
                }
            }
            for(int i = 0; i < 2; i++){
                for(int j = 0; j< tv.nj(); j++){
                    tv.endPoints()[i][j] = (eps) * tv.endPoints()[i][j] + (1.0-eps) * tvold.endPoints()[i][j];
                    if ( j  < (tv.nj()-1) ) tv.gamma()[i][j] = (eps) * tv.gamma()[i][j] + (1.0-eps) * tvold.gamma()[i][j];
                }
            }
        }
    }
    matchFirstWakePointsToSurface();
    for (int h = 0; h < (int) surfaces_.size(); h++){
        delete( olderSurfaces_[h] );
        olderSurfaces_[h] = oldSurfaces_[h];
        oldSurfaces_[h] =  new LiftingSurface( originalSurfaces[h] );
    }
    //Do time step to get velocities for PCC
}
void SimulationManager::solve(){
    if ( needsSolve_ ){
        int maxN = nOffset_.back();
        double* a = (double*) malloc( maxN*maxN*sizeof(double) );
        double* b = (double*) malloc( maxN*     sizeof(double) );
        memset(a, 0, maxN*maxN*sizeof(double));
        memset(b, 0, maxN*sizeof(double));

        fillRHS( b ); // Fill right hand side (b) of the linear system ( Ax = b )
        fillLHS( a ); // Fill left hand side (A) of the linear system  ( Ax = b )

        int n = maxN, nrhs = 1, lda = n, ldb = n, info;
        int* ipiv = (int*) malloc( maxN*sizeof(int) );
       
        //printMatrix( "LHS Matrix A", n, n,    a, lda );
        //printMatrix( "RHS Matrix B", n, nrhs, b, ldb );
       
        dgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );

        /* Check for the exact singularity */
        if( info > 0 ) {
                printf( "The diagonal element of the triangular factor of A,\n" );
                printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
                printf( "the solution could not be computed.\n" );
                exit( 1 );
        }
        //printMatrix( "Solution", n, nrhs, b, ldb );
        //printMatrix( "Details of LU factorization", n, n, a, lda );
        //print_int_vector( "Pivot indices", n, ipiv );
        
        //Loop back through to update gamma
        for( int h = 0; h < (int) surfaces_.size(); h++ ){
            LiftingSurface* s = surfaces_[h];
            HorseshoeLattice& hl = s->getHorseshoeLattice();
            std::fill(lastGamma_[h].begin(), lastGamma_[h].end(), 0.0 );
            std::fill(thisGamma_[h].begin(), thisGamma_[h].end(), 0.0 );
            //memset(&(lastGamma_[h][0]), lastGamma_[h].size()*sizeof( double ), 0 );
            //memset(&(thisGamma_[h][0]), thisGamma_[h].size()*sizeof( double ), 0 );
            for(int i = 0; i < s->nSpan(); i++){
                for(int j = 0; j < s->nChord(); j++){
                    //For every control point from every lifting surface
                    int n = hijToN(h,i,j);
                    lastGamma_[h][i] += hl.getGamma()[i][j]; //Accumulate last gammas for setting free wake BC
                    hl.getGamma()[i][j] = b[ipiv[n]-1]; 
                    thisGamma_[h][i] += hl.getGamma()[i][j]; //Accumulate these gammas for setting free wake BC
                }
            }
        }
        //Find net forces and moments on all lifting surfaces
        integrateForceAndMoment();
        free(a); free(b); free(ipiv);
    }
}
 
void SimulationManager::fillRHS( double* b ){
        //Outer loop indices and values have i appended, inner loop indices and values have j appended  
        for( int hi = 0; hi < (int) surfaces_.size(); hi++ ){
            LiftingSurface* si = surfaces_[hi];
            HorseshoeLattice& hli = si->getHorseshoeLattice();
            #pragma omp parallel for
            for(int ii = 0; ii < si->nSpan(); ii++){
                for(int ji = 0; ji < si->nChord(); ji++){
                    //For every control point from every lifting surface 
                    int ni = hijToN(hi,ii,ji);
                    Vec3D cpi  = hli.getControlPoints()[ii][ji]; //Save the location of this control point
                    Vec3D cpni  = hli.getControlPointNormals()[ii][ji]; //Save the location of this control point
                    Vec3D vInf = vInfinity( cpi );
                    b[ni] = -vInf.dot( cpni ); //Find the RHS (vInfinity dot control point normal)
                    
                    //Iterate through all attached vortex latices
                    for( int hj = 0; hj < (int)surfaces_.size(); hj++ ){
                        LiftingSurface* sj = surfaces_[hj];
                        if ( sj->freeWake() ){
                            VortexLattice& vlj = sj->getVortexLattice();
                            b[ni] -= vlj.calcInducedVelocity( cpi, 1 ).dot( cpni );
                        }
                        if ( sj->freeTipVortex() ){
                            TipFilament& tf = sj->getTipFilament();
                            b[ni] -= tf.calcInducedVelocity( cpi ).dot( cpni );
                        }
                    }//hj
                }//ji
            }//ii
        }//hi
}

void SimulationManager::fillLHS( double* a ){
        //Outer loop indices and values have i appended, inner loop indices and values have j appended  
        for( int hi = 0; hi < (int) surfaces_.size(); hi++ ){
            LiftingSurface* si = surfaces_[hi];
            HorseshoeLattice& hli = si->getHorseshoeLattice();
            #pragma omp parallel for
            for(int ii = 0; ii < si->nSpan(); ii++){
                for(int ji = 0; ji < si->nChord(); ji++){
                    //For every control point from every lifting surface
                    int ni = hijToN(hi,ii,ji);
                    Vec3D cpi  = hli.getControlPoints()[ii][ji]; //Save the location of this control point
                    Vec3D cpni  = hli.getControlPointNormals()[ii][ji]; //Save the location of this control point
                    //Iterate through all filaments 
                    for( int hj = 0; hj < (int)surfaces_.size(); hj++ ){
                        LiftingSurface* sj = surfaces_[hj];
                        HorseshoeLattice& hlj = sj->getHorseshoeLattice();
                        for(int ij = 0; ij < sj->nSpan(); ij++){
                            for(int jj = 0; jj < sj->nChord(); jj++){
                                //For every horseshoe segment from every lifting surface
                                int nj = hijToN(hj,ij,jj);
                                int superN = nInJToSuperN( ni, nj );
                                //Find the influence coefficient from vortex (hj,ij,jj) to control point (hi,ii,ji) 
                                a[superN] += cpni.dot( hlj.calcInfluenceCoefficient( cpi, hlj.ijToN(ij,jj) ) ); 
                                //printf("set a at %i, %i [%i] = %16.16f \n", ni, nj, superN, a[superN]);
                            }//jj
                        }//ij
                    }//hj
                }//ji
            }//ii
        }//hi
}

double SimulationManager::netLiftDirect(){
    double netLift;
    for( int h = 0; h < (int) surfaces_.size(); h++ ){
        LiftingSurface* s = surfaces_[h];
        HorseshoeLattice& hl = s->getHorseshoeLattice();
        for(int i = 0; i < s->nSpan(); i++){
            for(int j = 0; j < s->nChord(); j++){
                //For every control point from every lifting surface
                Vec3D vInf  = vInfinity( hl.getControlPoints()[i][j] );
                netLift += hl.getGamma()[i][j]*vInf.magnitude()*hl.dSpan(i,j);
            }
        }
    }
                    
   return netLift; 
}

void SimulationManager::integrateForceAndMoment(){
    Vec3D netForce = Vec3D();
    Vec3D netMoment = Vec3D();
    for( int h = 0; h < (int) surfaces_.size(); h++ ){
        LiftingSurface* s = surfaces_[h];
        HorseshoeLattice& hl = s->getHorseshoeLattice();
        for(int i = 0; i < s->nSpan(); i++){
            Vec3D chordwiseSum = Vec3D();
            for(int j = 0; j < s->nChord(); j++){
                Vec3D vInduced = Vec3D(0.0, 0.0, 0.0);
                Vec3D thisGamma = hl.gammaVector(i,j);
                Vec3D thisGammaCenter = hl.gammaCenterPoint(i,j);
                Vec3D vInf  = vInfinity( thisGammaCenter );
                //Iterate through all filaments 
                for( int hj = 0; hj < (int) surfaces_.size(); hj++ ){
                    Vec3D vIndH = Vec3D(0.0, 0.0, 0.0);
                    LiftingSurface* sj = surfaces_[hj];
                    //HorseshoeLattice& hlj = sj->getHorseshoeLattice();
                    vIndH = sj->calcInducedVelocity( thisGammaCenter );
                    //vInduced += hlj.calcInducedVelocity( thisGammaCenter );
                    vInduced += vIndH;
                }//hj
                Vec3D force;
                if ( refSurf_.pgCorrection ){
                    double beta = 1.0/sqrt(1.0 - (vInf+vInduced).magnitude()*(vInf+vInduced).magnitude()/(refSurf_.vMach * refSurf_.vMach));
                    force = (vInf+vInduced).cross( thisGamma )/beta;
                } else {
                    force = (vInf+vInduced).cross( thisGamma );
                }
                netForce += force;
                netMoment += force.cross(thisGammaCenter);
                chordwiseSum += force;
            }
            s->spanwiseForce()[i] = chordwiseSum/(hl.dSpan(i,0));
        }
    }
    fomo_.bodyForce = netForce;
    fomo_.bodyMoment = netMoment;
    fomo_.bodyForceCoeff = netForce/(1.0/2.0 * refV_ * refV_ * refSurf_.S );
    fomo_.bodyMomentCoeff = netMoment/(1.0/2.0 * refV_ * refV_ * refSurf_.S * refSurf_.b );
    Vec3D uX = Vec3D(1.0, 0.0, 0.0);
    Vec3D uV = globalLinearVelocity_.norm();
    Vec3D uL = uV.rotate( Vec3D(), Vec3D(0.0, 1.0, 0.0), -M_PI / 2.0 );
    Vec3D binormal = uV.cross( uL );
    fomo_.aeroForce = Vec3D( fomo_.bodyForce.dot(uV), fomo_.bodyForce.dot(binormal), fomo_.bodyForce.dot(uL) );
    fomo_.aeroForceCoeff = fomo_.aeroForce /  (1.0/2.0 * refV_ * refV_ * refSurf_.S );
}

void SimulationManager::calculateWakeVelocities(){
    //For every surface
    for( int h = 0; h < (int) surfaces_.size(); h++ ){
        LiftingSurface* s = surfaces_[h];
        if  ( s->freeWake() ){ //If surface is a free wake case we need to calculate wake velocity
            VortexLattice& vl = s->getVortexLattice();
            //For end points along span and wake
            for(int i = 0; i < s->nSpan()+1; i++){ 
                for(int j = 0; j < s->nWake(); j++){ 
                    Vec3D thisEndPoint = vl.endPoints()[i][j];
                    Vec3D vInduced = Vec3D();
                    //Iterate through all other surfaces to find net induced velocity 
                    for( int hj = 0; hj < (int) surfaces_.size(); hj++ ){
                        LiftingSurface* sj = surfaces_[hj];
                        if (j > 0) vInduced += sj->calcInducedVelocity( thisEndPoint );
                    }//hj
                    Vec3D vInf  = vInfinity( thisEndPoint );
                    //Vec3D vInf  = vInfinityLinearOnly();
                    vl.endPointVelocity()[i][j] = vInduced + vInf;
                }
            }
        }
        if  ( s->freeTipVortex() ){ //If surface has a free tip filament we need to calculate wake velocity
            TipFilament& tf = s->getTipFilament();
            //For two tip vortex filaments
            for(int i = 0; i < 2; i++){ 
                for(int j = 0; j < s->nFilament(); j++){ 
                    Vec3D thisEndPoint = tf.endPoints()[i][j];
                    Vec3D vInduced = Vec3D();
                    //Iterate through all other surfaces to find net induced velocity 
                    for( int hj = 0; hj < (int) surfaces_.size(); hj++ ){
                        LiftingSurface* sj = surfaces_[hj];
                        if (j > 0) vInduced += sj->calcInducedVelocity( thisEndPoint );
                    }//hj
                    Vec3D vInf  = vInfinity( thisEndPoint );
                    //Vec3D vInf  = vInfinityLinearOnly();
                    tf.endPointVelocity()[i][j] = vInduced + vInf;
                }
            }
        }
    }
}

void SimulationManager::advectWake(){
    //For every surface
    for( int h = 0; h < (int) surfaces_.size(); h++ ){
        LiftingSurface* s = surfaces_[h];
        if  ( s->freeWake() ){ //If surface is a free wake case we need to calculate wake velocity
            VortexLattice& vl = s->getVortexLattice();
            vl.advectAndRotate( dt_, globalRotationAxis_, globalRotationRate_ );
            TipFilament& tf = s->getTipFilament();
            tf.fixToWake( vl );
            tf.advectAndRotate( dt_, globalRotationAxis_, globalRotationRate_ );
        }
    }
}

void SimulationManager::advectWakePC2B(){
    //For every surface
    for( int h = 0; h < (int) surfaces_.size(); h++ ){
        LiftingSurface* s = surfaces_[h];
        if  ( s->freeWake() ){ //If surface is a free wake case we need to calculate wake velocity
            VortexLattice& vl = s->getVortexLattice();
            vl.advectPC2B( dt_, globalRotationAxis_, globalRotationRate_, oldSurfaces_[h]->getVortexLattice(), olderSurfaces_[h]->getVortexLattice() );
            TipFilament& tf = s->getTipFilament();
            tf.fixToWake( vl );
            tf.advectPC2B( dt_, globalRotationAxis_, globalRotationRate_, oldSurfaces_[h]->getTipFilament(), olderSurfaces_[h]->getTipFilament() );
        }
    }
}

void SimulationManager::matchFirstWakePointsToSurface(){
    for (int h = 0; h < (int) surfaces_.size(); h++){
        LiftingSurface* ls = surfaces_[h];
        if (ls->freeWake()){
            ls->getVortexLattice().fixToTrailingEdge( ls->getHorseshoeLattice() );
        }
    }
}

void SimulationManager::fillWakeBC(){
    //For every surface
    for( int h = 0; h < (int) surfaces_.size(); h++ ){
        LiftingSurface* s = surfaces_[h];
        if  ( s->freeWake() ){ //If surface is a free wake case we need to calculate wake velocity
            VortexLattice& vl = s->getVortexLattice();
            vl.fixToTrailingEdge( s->getHorseshoeLattice() );
            //For end points along span and wake
            for(int i = 0; i < s->nSpan()+1; i++){
                if ( i < s->nSpan() ) vl.gammaI()[i][0] = (lastGamma_[h][i] - thisGamma_[h][i])*.25;
                if ( i > 0 && i < s->nSpan() ) {   
                    vl.gammaJ()[i][0] = lastGamma_[h][i-1] - lastGamma_[h][i]; 
                } else if ( i == s->nSpan() ) {   
                    vl.gammaJ()[i][0] = lastGamma_[h][i-1]; 
                } else if ( i == 0 ) {
                    vl.gammaJ()[i][0] = -lastGamma_[h][i];
                }
            }
            TipFilament& tf = s->getTipFilament();
            tf.fixToWake( vl );
        }
    }
}

double SimulationManager::netDrag(){
    return fomo_.aeroForce.x; 
}

double SimulationManager::netLift(){
    return fomo_.aeroForce.z;
}

void SimulationManager::setReferenceSurface( ReferenceSurface s ){
    refSurf_ = s;
}

void SimulationManager::setReferenceVelocity( double v ){
    refV_ = v;
}

void SimulationManager::setDt( double dt ){
    dt_ = dt;
}

ReferenceSurface SimulationManager::referenceSurface(){
    return refSurf_;
}

double SimulationManager::referenceVelocity(){
    return refV_;
}

double SimulationManager::dt(){
    return dt_;
}

ForcesAndMoments SimulationManager::forcesAndMoments(){
    return fomo_;
}

LiftingSurface& SimulationManager::getSurface(int i){
    return *(surfaces_[i]);
}

int SimulationManager::getNSurfaces(){
    return surfaces_.size();
}

Vec3D SimulationManager::vInfinity( Vec3D p ){
    Vec3D vInf = globalLinearVelocity_; 
    //Vec3d vRot = globablRotationRate_ * globalRotationAxis_.cross(p);
    if( fabs(globalRotationRate_) > 0.0 ){
        Vec3D vRot1 = globalRotationRate_ * globalRotationAxis_.cross( p );
        vInf += vRot1;
    }
    return vInf;
}

Vec3D SimulationManager::vInfinityLinearOnly(){
    Vec3D vInf = globalLinearVelocity_; 
    return vInf;
}

int SimulationManager::hijToN( int h, int i, int j ){
    int n = 0;
    n += nOffset_[h];
    n += surfaces_[h]->getHorseshoeLattice().ijToN(i,j);
    return n; 
}

int SimulationManager::nInJToSuperN( int ni, int nj ){
    int maxN = nOffset_.back();
    int superN = nj*maxN+ni; 
    return superN;
}

std::tuple<int, int, int> SimulationManager::hijFromN( int n ){
    int h = 0;
    while ( nOffset_[h] < n ) h++;
    n -= nOffset_[h];
    std::pair<int, int> ij = surfaces_[h]->getHorseshoeLattice().ijFromN(n);
    return std::make_tuple(n, ij.first, ij.second); 
}

void SimulationManager::printState(){
    double Cl =  fomo_.aeroForceCoeff.z;
    double Cd =  fomo_.aeroForceCoeff.x;
    double CdIdeal = Cl*Cl/ (M_PI * refSurf_.S );
    printf("The net lift was found to be %8.8f\n", netLift() );
    printf("The net lift was found to be %8.8f  by direct integration\n", netLiftDirect() );
    printf("Cl was found to be %8.8f\n", Cl );
    printf("The net drag was found to be %8.8f\n", netDrag() );
    printf("Cd was found to be %8.8f\n", Cd );
    printf("Cd ideal is %8.8f\n", CdIdeal );
    printf("Oswald efficient factor of %8.8f\n", CdIdeal/Cd );
    printf("forces\n");
    fomo_.bodyForce.printState();
    fomo_.aeroForce.printState();
    fomo_.aeroForceCoeff.printState();
}

void printMatrix( char const* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %17.17f", a[i*lda+j] );
                printf( "\n" );
        }
}

void print_int_vector( char const* desc, int n, int* a ) {
        int j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
        printf( "\n" );
}

