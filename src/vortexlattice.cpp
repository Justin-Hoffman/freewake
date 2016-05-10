#include "vortexlattice.h"
#include "vortexmath.h"
#include "stdio.h"

VortexLattice::VortexLattice() : ni_( 2 ), nj_( 2 ), 
                                 endPoints_( 2, std::vector<Vec3D>( 2, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 endPointV_( 2, std::vector<Vec3D>( 2, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 gammaI_( 1, std::vector<double>( 2, 0.0  ) ) , 
                                 rcI_( 1, std::vector<double>( 2, 1.E-2  ) ) , 
                                 gammaJ_( 2, std::vector<double>( 1, 0.0  ) ) ,
                                 rcJ_( 2, std::vector<double>( 1, 1.E-1 ) ) {
 
}

VortexLattice::VortexLattice( int ni, int nj ) : 
                                 ni_( ni ), nj_( nj ),
                                 endPoints_( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 endPointV_( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 gammaI_( ni-1, std::vector<double>( nj, 0.0  ) ), 
                                 rcI_( ni-1, std::vector<double>( nj, 1.E-2  ) ), 
                                 gammaJ_( ni, std::vector<double>( nj-1, 0.0  ) ),
                                 rcJ_( ni, std::vector<double>( nj-1, 1.E-1  ) ){
 
}

VortexLattice::VortexLattice( const VortexLattice &v ) : 
                                 ni_( v.ni_ ), nj_( v.nj_ ),
                                 endPoints_( v.endPoints_ ),
                                 endPointV_( v.endPointV_ ),
                                 gammaI_( v.gammaI_ ), 
                                 rcI_( v.rcI_ ), 
                                 gammaJ_( v.gammaJ_ ),
                                 rcJ_( v.rcJ_ ){
 
}

VortexLattice::~VortexLattice() {
    
}

Vec3D VortexLattice::calcInfluenceCoefficient( Vec3D p, int n ){
    Vec3D vOut = Vec3D(p);
    std::pair<int, int> ij = ijFromN( n );
    int i = ij.first; 
    int j = ij.second;
    Vec3D r1, r2, aOut;
    aOut = Vec3D(0.0, 0.0, 0.0);
    if ( i < ni_ ){ 
        r1 = (endPoints_[i][j]);
        r2 = (endPoints_[i+1][j]);
        aOut = BiotSavart(r1, r2, rcI_[i][j] );
    }
    if ( j < nj_ ){
        r1 = (endPoints_[i][j]);
        r2 = (endPoints_[i][j+1]);
        aOut += BiotSavart(r1, r2, rcJ_[i][j] );
    } 

    return aOut;
}

Vec3D VortexLattice::calcInducedVelocity( Vec3D p, int jStart ){
    Vec3D r1, r2, aOut = Vec3D();
    double vx = 0.0, vy = 0.0, vz = 0.0;
    #pragma omp parallel for private( r1, r2) reduction(+:vx,vy,vz)
    for ( int i = 0; i < ni_; i++){
        Vec3D aTmp = Vec3D(0.0, 0.0, 0.0);//Hackery because I can't use openMP reduce on a class
        for (int j = jStart; j < nj_; j++){
            if (i < ni_-1) {
                //Contribution of spanwise (i) edge filaments
                r1 = (endPoints_[i  ][j] - p);
                r2 = (endPoints_[i+1][j] - p);
                aTmp += BiotSavart( r1, r2, rcI_[i][j] )*gammaI_[i][j];
            }
            if (j < nj_-1) {
                //Contribution of chordwise (j) edge filaments
                r1 = (endPoints_[i][j  ] - p);
                r2 = (endPoints_[i][j+1] - p);
                aTmp += BiotSavart( r1, r2, rcJ_[i][j] )*gammaJ_[i][j];
            }    
        }
        vx += aTmp.x;
        vy += aTmp.y;
        vz += aTmp.z;
    }
    aOut.x = vx; aOut.y = vy; aOut.z = vz;
    return aOut;
}

void VortexLattice::advect( double dt ){
    for ( int i = ni_-1; i > -1; i--){
        for (int j = nj_-1; j > 0; j--){
            endPoints_[i][j] = endPoints_[i][j-1] + (endPointV_[i][j-1])*dt;
            if ( i < ni_-1 ){
                gammaI_[i][j] = gammaI_[i][j-1];
                rcI_[i][j] = VortexCoreGrowth( rcI_[i][j-1], dt );
            }
            if ( j < nj_-1 ){
                gammaJ_[i][j] = gammaJ_[i][j-1];
                rcJ_[i][j] = VortexCoreGrowth( rcJ_[i][j-1], dt );
            }
        }
    }    
}

void VortexLattice::advectAndRotate( double dt, Vec3D axis, double omega ){
    for ( int i = ni_-1; i > -1; i--){
        for (int j = nj_-1; j > 0; j--){
            endPoints_[i][j] = endPoints_[i][j-1].rotate(Vec3D(0.0, 0.0, 0.0), axis, omega*dt) + (endPointV_[i][j-1])/1.0*dt;
            
            if ( i < ni_-1 ){
                gammaI_[i][j] = gammaI_[i][j-1];
                rcI_[i][j] = VortexCoreGrowth( rcI_[i][j-1], dt );
            }
            if ( j < nj_-1 ){
                gammaJ_[i][j] = gammaJ_[i][j-1];
                rcJ_[i][j] = VortexCoreGrowth( rcJ_[i][j-1], dt );
            }
        }
    }    
}

void VortexLattice::advectPCC( double dt, Vec3D axis, double omega ){
    for ( int i = ni_-1; i > -1; i--){
        Vec3D pjp1tp0 = Vec3D(endPoints_[i][0]), pjp0tp0 = Vec3D(endPoints_[i][0]); //Points at this time at 2 wake ages
        for (int j = 0; j < nj_-1; j++){
            pjp0tp0 = Vec3D(pjp1tp0);     //Point at j for current time is whatever we had saved
            pjp1tp0 = endPoints_[i][j+1]; //Save point at j+1 for current time
            endPoints_[i][j+1] = ( (endPointV_[i][j] + endPointV_[i][j+1]) / 1.0 * dt - ( 1.0-1.0) *       endPoints_[i][j] // j+1 t+1 is RHS
                                                                                      - (-1.0-1.0) *   pjp0tp0              - (-1.0+1.0) *   pjp1tp0 ) / (1.0+1.0);
        }
    }    
    for ( int i = ni_-1; i > -1; i--){
        for (int j = nj_-1; j > 0; j--){
            if ( i < ni_-1 ){
                gammaI_[i][j] = gammaI_[i][j-1];
                rcI_[i][j] = VortexCoreGrowth( rcI_[i][j-1], dt );
            }
            if ( j < nj_-1 ){
                gammaJ_[i][j] = gammaJ_[i][j-1];
                rcJ_[i][j] = VortexCoreGrowth( rcJ_[i][j-1], dt );
            }
        }
    }
}

void VortexLattice::advectPC2B( double dt, Vec3D axis, double omega, VortexLattice& old, VortexLattice& older ){
    for ( int i = ni_-1; i > -1; i--){
        Vec3D pjp1tp0 = Vec3D(endPoints_[i][0]), pjp0tp0 = Vec3D(endPoints_[i][0]); //Points at this time at 2 wake ages
        for (int j = 0; j < nj_-1; j++){
            pjp0tp0 = Vec3D(pjp1tp0);     //Point at j for current time is whatever we had save
            pjp1tp0 = endPoints_[i][j+1]; //Save point at j+1 for current time
            endPoints_[i][j+1] = ( (endPointV_[i][j] + endPointV_[i][j+1])/1.0 * dt - ( 3.0/4.0-1.0) *       endPoints_[i][j] // j+1 t+1 is RHS 
                                                                                    - (-1.0/4.0-1.0) *   pjp0tp0              - (-1.0/4.0+1.0) *   pjp1tp0
                                                                                    - (-3.0/4.0+0.0) *   old.endPoints_[i][j] - (-3.0/4.0+0.0) *   old.endPoints_[i][j+1]
                                                                                    - ( 1.0/4.0+0.0) * older.endPoints_[i][j] - ( 1.0/4.0+0.0) * older.endPoints_[i][j+1]) / ( 3.0/4.0+1.0 );
        }
    }
    for ( int i = ni_-1; i > -1; i--){
        for (int j = nj_-1; j > 0; j--){
            if ( i < ni_-1 ){
                gammaI_[i][j] = gammaI_[i][j-1];
                rcI_[i][j] = VortexCoreGrowth( rcI_[i][j-1], dt );
            }
            if ( j < nj_-1 ){
                gammaJ_[i][j] = gammaJ_[i][j-1];
                rcJ_[i][j] = VortexCoreGrowth( rcJ_[i][j-1], dt );
            }
        }
    }    
}

void VortexLattice::fixToTrailingEdge( HorseshoeLattice &h ){
    int maxj = h.nj();
    for ( int i = 0; i < ni_; i++){
        endPoints_[i][0] = h.getEndPoints()[i][maxj];
        h.getTrailedPoints()[i] = Vec3D(endPoints_[i][1]);
    }
}

int VortexLattice::ni(){
    return ni_;
}

int VortexLattice::nj(){
    return nj_;
}

std::vector<std::vector<Vec3D>>& VortexLattice::endPoints(){
    return endPoints_;
} 

std::vector<std::vector<Vec3D>>& VortexLattice::endPointVelocity(){
    return endPointV_;
} 

std::vector<std::vector<double>>& VortexLattice::gammaI(){
    return gammaI_;
} 

std::vector<std::vector<double>>& VortexLattice::gammaJ(){
    return gammaJ_;
} 

std::vector<std::vector<double>>& VortexLattice::rcI(){
    return rcI_;
} 

std::vector<std::vector<double>>& VortexLattice::rcJ(){
    return rcJ_;
} 


std::pair<int, int> VortexLattice::ijFromN( int n ){
    std::pair<int, int> ij = std::pair<int, int>();
    ij.first = n/ni_;
    ij.second = n%(ni_);

    return ij;
}

void VortexLattice::initializeToHelix( Vec3D axis, double dTheta, double dZ ){
    for ( int i = 0; i < ni_; i++){
        for (int j = 1; j < nj_; j++){
            endPoints_[i][j] = endPoints_[i][j-1].rotate( Vec3D(), axis, dTheta ) - axis.norm()*dZ;
        }
    }    
}

void VortexLattice::printState(){
    for (int i = 0; i < ni_; i++){
        for(int j = 0; j < nj_; j++){
            printf("loop at i = %i\n endPoint at",i);   
            endPoints_[i][j].printState();      
            endPointV_[i][j].printState();     
    
        }
    }
    for (int i = 0; i < ni_-1; i++){
        for(int j = 0; j < nj_-1; j++){
            printf("loop at i = %i rcI/J = %5.5e / %5.5e \n",i, rcI_[i][j], rcJ_[i][j]); 
        }
    }  
}
