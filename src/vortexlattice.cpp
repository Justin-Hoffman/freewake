#include "vortexlattice.h"
#include "vortexmath.h"

VortexLattice::VortexLattice() : ni_( 2 ), nj_( 2 ), 
                                 endPoints( 2, std::vector<Vec3D>( 2, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 gammaI( 2, std::vector<double>( 2, 0.0  ) ) , 
                                 gammaJ( 2, std::vector<double>( 2, 0.0  ) ) {
 
}

VortexLattice::VortexLattice( int ni, int nj ) : 
                                 ni_( ni ), nj_( nj ),
                                 endPoints( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 gammaI( ni-1, std::vector<double>( nj, 0.0  ) ), 
                                 gammaJ( ni, std::vector<double>( nj-1, 0.0  ) ){
 
}

VortexLattice::VortexLattice( const VortexLattice &v ) : 
                                 ni_( v.ni_ ), nj_( v.nj_ ),
                                 endPoints( v.endPoints ),
                                 gammaI( v.gammaI ), 
                                 gammaJ( v.gammaJ ){
 
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
        r1 = (endPoints[i][j]);
        r2 = (endPoints[i+1][j]);
        aOut = BiotSavart(r1, r2, 0.0);
    }
    if ( j < nj_ ){
        r1 = (endPoints[i][j]);
        r2 = (endPoints[i][j+1]);
        aOut += BiotSavart(r1, r2, 0.0);
    } 

    return aOut;
}

Vec3D VortexLattice::calcInducedVelocity( Vec3D p){
    Vec3D vOut = Vec3D(p);

    return vOut;
}

int VortexLattice::ni(){
    return ni_;
}

int VortexLattice::nj(){
    return nj_;
}

std::vector<std::vector<Vec3D>>& VortexLattice::getEndPoints(){
    return endPoints;
} 

std::vector<std::vector<double>>& VortexLattice::getGammaI(){
    return gammaI;
} 

std::vector<std::vector<double>>& VortexLattice::getGammaJ(){
    return gammaJ;
} 

std::pair<int, int> VortexLattice::ijFromN( int n ){
    std::pair<int, int> ij = std::pair<int, int>();
    ij.first = n/ni_;
    ij.second = n%(ni_);

    return ij;
}
