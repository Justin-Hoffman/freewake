#include <cmath>
#include "simulationmanager.h"
#include "lapacke.h"
#include "string.h"

SimulationManager::SimulationManager() : needsSolve_( true ), surfaces_(), nOffset_(1,0), a(NULL), b(NULL), refSurf_( ReferenceSurface(1.0, 1.0, 1.0) ), refV_(1.0), 
                                         globalLinearVelocity_(), globalRotationAxis_(0.0,0.0,1.0), globalRotationRate_(0.0), 
                                         fomo_() {

}

SimulationManager::SimulationManager( SimulationManager& s ) : needsSolve_( s.needsSolve_ ),  surfaces_(s.surfaces_), nOffset_(s.nOffset_), a(NULL), b(NULL), refSurf_( s.refSurf_ ), refV_( s.refV_ ), 
                                                               globalLinearVelocity_( s.globalLinearVelocity_), globalRotationAxis_( s.globalRotationAxis_ ), globalRotationRate_( s.globalRotationRate_),
                                                               fomo_( s.fomo_ ) {

}

SimulationManager::~SimulationManager(){
}

void SimulationManager::addSurface( LiftingSurface* s){
    surfaces_.push_back( s );
    nOffset_.push_back( s->getLattice().maxN() + nOffset_.back() );
}

void SimulationManager::setGlobalLinearVelocity( Vec3D v ){
    globalLinearVelocity_ = v;
}

Vec3D SimulationManager::getGlobalLinearVelocity(){
    return globalLinearVelocity_;
}

extern void printMatrix( char const* desc, int m, int n, double* a, int lda );
extern void print_int_vector( char const* desc, int n, int* a );
 
void SimulationManager::solve(){
    if ( needsSolve_ ){
        int maxN = nOffset_.back();
        double* a = (double*) malloc( maxN*maxN*sizeof(double) );
        double* b = (double*) malloc( maxN*     sizeof(double) );
        memset(a, 0, maxN*maxN*sizeof(double));
      
        //Outer loop indices and values have i appended, inner loop indices and values have j appended  
        for( int hi = 0; hi < (int) surfaces_.size(); hi++ ){
            LiftingSurface* si = surfaces_[hi];
            HorseshoeLattice& hli = si->getLattice();
            #pragma omp parallel for
            for(int ii = 0; ii < si->nSpan(); ii++){
                for(int ji = 0; ji < si->nChord(); ji++){
                    //For every control point from every lifting surface
            
                    int ni = hijToN(hi,ii,ji);
                    Vec3D cpi  = hli.getControlPoints()[ii][ji]; //Save the location of this control point
                    Vec3D cpni  = hli.getControlPointNormals()[ii][ji]; //Save the location of this control point
                    Vec3D vInf = vInfinity( cpi );
                    b[ni] = -vInf.dot( cpni ); //Find the RHS (vInfinity dot control point normal)
                    //printf("set b at %i = %16.16f \n", ni, b[ni]);
                       
                   
                    //Iterate through all filaments 
                    for( int hj = 0; hj < (int)surfaces_.size(); hj++ ){
                        LiftingSurface* sj = surfaces_[hj];
                        HorseshoeLattice& hlj = sj->getLattice();
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

        /* Locals */
        int n = maxN, nrhs = 1, lda = n, ldb = nrhs, info;
        /* Local arrays */
        int* ipiv = (int*) malloc( maxN*sizeof(int) );
       
        //printMatrix( "LHS Matrix A", n, n,    a, lda );
        //printMatrix( "RHS Matrix B", n, nrhs, b, ldb );
       
        /* Solve the equations A*X = B */
        info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );

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
            HorseshoeLattice& hl = s->getLattice();
            for(int i = 0; i < s->nSpan(); i++){
                for(int j = 0; j < s->nChord(); j++){
                    //For every control point from every lifting surface
                    int n = hijToN(h,i,j);
                    hl.getGamma()[i][j] = b[ipiv[n]-1]; 
                }
            }
        }
        integrateForceAndMoment();
        free(a); free(b);
    }
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

double SimulationManager::netLiftDirect(){
    double netLift;
    for( int h = 0; h < (int) surfaces_.size(); h++ ){
        LiftingSurface* s = surfaces_[h];
        HorseshoeLattice& hl = s->getLattice();
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
        HorseshoeLattice& hl = s->getLattice();
        for(int i = 0; i < s->nSpan(); i++){
            for(int j = 0; j < s->nChord(); j++){
                Vec3D vInduced = Vec3D(0.0, 0.0, 0.0);
                Vec3D thisGamma = hl.gammaVector(i,j);
                Vec3D thisGammaCenter = hl.gammaCenterPoint(i,j);
                Vec3D vInf  = vInfinity( thisGammaCenter );
                //Iterate through all filaments 
                for( int hj = 0; hj < (int) surfaces_.size(); hj++ ){
                    Vec3D vIndH = Vec3D(0.0, 0.0, 0.0);
                    LiftingSurface* sj = surfaces_[hj];
                    HorseshoeLattice& hlj = sj->getLattice();
                    vIndH = hlj.calcInducedVelocity( thisGammaCenter );
                    //vInduced += hlj.calcInducedVelocity( thisGammaCenter );
                    vInduced += vIndH;
                }//hj
                Vec3D force = (vInf+vInduced).cross( thisGamma );
                netForce += force;
                netMoment += force.cross(thisGammaCenter);
            }
        }
    }
    fomo_.bodyForce = netForce;
    fomo_.bodyMoment = netMoment;
    fomo_.bodyForceCoeff = netForce/(1.0/2.0 * refV_ * refV_ * refSurf_.S );
    fomo_.bodyMomentCoeff = netForce/(1.0/2.0 * refV_ * refV_ * refSurf_.S );
    Vec3D uX = Vec3D(1.0, 0.0, 0.0);
    Vec3D uV = globalLinearVelocity_.norm();
    Vec3D uL = uV.rotate( Vec3D(), Vec3D(0.0, 1.0, 0.0), -M_PI / 2.0 );
    Vec3D binormal = uV.cross( uL );
    fomo_.aeroForce = Vec3D( fomo_.bodyForce.dot(uV), fomo_.bodyForce.dot(binormal), fomo_.bodyForce.dot(uL) );
    fomo_.aeroForceCoeff = fomo_.aeroForce /  (1.0/2.0 * refV_ * refV_ * refSurf_.S );
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

ReferenceSurface SimulationManager::referenceSurface(){
    return refSurf_;
}

double SimulationManager::referenceVelocity(){
    return refV_;
}

Vec3D SimulationManager::vInfinity( Vec3D p ){
    Vec3D vInf = globalLinearVelocity_;  
    vInf += globalRotationRate_ * globalRotationAxis_.cross(p);
    return vInf;
}

int SimulationManager::hijToN( int h, int i, int j ){
    int n = 0;
    n += nOffset_[h];
    n += surfaces_[h]->getLattice().ijToN(i,j);
    return n; 
}

int SimulationManager::nInJToSuperN( int ni, int nj ){
    int maxN = nOffset_.back();
    int superN = ni*maxN+nj; 
    return superN;
}

std::tuple<int, int, int> SimulationManager::hijFromN( int n ){
    int h = 0;
    while ( nOffset_[h] < n ) h++;
    n -= nOffset_[h];
    std::pair<int, int> ij = surfaces_[h]->getLattice().ijFromN(n);
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
