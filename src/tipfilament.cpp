#include <cmath>
#include "tipfilament.h"
#include "vortexmath.h"
#include "stdio.h"

TipFilament::TipFilament() : ni_( 2 ), nj_( 2 ), 
                                 endPoints_( 2, std::vector<Vec3D>( 2, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 endPointV_( 2, std::vector<Vec3D>( 2, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 gamma_( 2, std::vector<double>( 1, 0.0  ) ) ,
                                 rc_( 2, std::vector<double>( 1, 5.E-1 ) ) {
 
}

TipFilament::TipFilament( int ni, int nj ) : 
                                 ni_( ni ), nj_( nj ),
                                 endPoints_( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 endPointV_( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 gamma_( ni, std::vector<double>( nj-1, 0.0  ) ), 
                                 rc_( ni, std::vector<double>( nj-1, 5.E-1  ) ){
 
}

TipFilament::TipFilament( const TipFilament &tf ) : 
                                 ni_( tf.ni_ ), nj_( tf.nj_ ),
                                 endPoints_( tf.endPoints_ ),
                                 endPointV_( tf.endPointV_ ),
                                 gamma_( tf.gamma_ ),
                                 rc_( tf.rc_ ){
}

TipFilament::~TipFilament() {
    
}

Vec3D TipFilament::calcInfluenceCoefficient( Vec3D p, int n ){
    Vec3D vOut = Vec3D(p);
    std::pair<int, int> ij = ijFromN( n );
    int i = ij.first; 
    int j = ij.second;
    Vec3D r1, r2, aOut;
    aOut = Vec3D(0.0, 0.0, 0.0);
    if ( j < nj_-1 ){
        r1 = (endPoints_[i][j]);
        r2 = (endPoints_[i][j+1]);
        aOut += BiotSavart(r1, r2, rc_[i][j] );
    } 

    return aOut;
}

Vec3D TipFilament::calcInducedVelocity( Vec3D p, int jStart ){
    Vec3D r1, r2, aOut = Vec3D();
    double vx = 0.0, vy = 0.0, vz = 0.0;
    #pragma omp parallel for private( r1, r2) reduction(+:vx,vy,vz)
    for ( int i = 0; i < ni_; i++){
        Vec3D aTmp = Vec3D(0.0, 0.0, 0.0);//Hackery because I can't use openMP reduce on a class
        for (int j = jStart; j < nj_; j++){
            if (j < nj_-1) {
                //Contribution of chordwise (j) edge filaments
                r1 = (endPoints_[i][j  ] - p);
                r2 = (endPoints_[i][j+1] - p);
                aTmp += BiotSavart( r1, r2, rc_[i][j] )*gamma_[i][j];
            }    
        }
        vx += aTmp.x;
        vy += aTmp.y;
        vz += aTmp.z;
    }
    aOut.x = vx; aOut.y = vy; aOut.z = vz;
    return aOut;
}

void TipFilament::advect( double dt ){
    for ( int i = ni_-1; i > -1; i--){
        for (int j = nj_-1; j > 0; j--){
            endPoints_[i][j] = endPoints_[i][j-1] + endPointV_[i][j-1]*dt;
            if ( j < nj_-1 ){
                gamma_[i][j] = gamma_[i][j-1];
                rc_[i][j] = VortexCoreGrowth( rc_[i][j-1], dt );
            }
        }
    }    
}

void TipFilament::advectAndRotate( double dt, Vec3D axis, double omega ){
    for ( int i = ni_-1; i > -1; i--){
        for (int j = nj_-1; j > 0; j--){
            endPoints_[i][j] = endPoints_[i][j-1].rotate(Vec3D(0.0, 0.0, 0.0), axis, omega*dt) + endPointV_[i][j-1]*dt;
            if ( j < nj_-1 ){
                gamma_[i][j] = gamma_[i][j-1];
                rc_[i][j] = VortexCoreGrowth( rc_[i][j-1], dt );
            }
        }
    }    
}

void TipFilament::fixToWake( VortexLattice &vl ){
    int maxi = vl.ni();
    int maxj = vl.nj();
    double nGamma = 0;
    double pGamma = 0;
    int maxGammaI = 0;
    for (int i = 0; i < maxi; i++){ 
        if ( vl.gammaJ()[i][maxj-2] < 0.0 ){
            nGamma += vl.gammaJ()[i][maxj-2];
            maxGammaI = i; 
        } else {
            pGamma += vl.gammaJ()[i][maxj-2];
        }
    }
    double den = 0.0, numx = 0.0, numy = 0.0, numz=0.0;
    double thisGamma = 0.0;
    for(int i = 0; i < maxGammaI/2; i++){
        thisGamma = 0.0;
        if ( i > 0 ) { thisGamma += vl.gammaJ()[i-1][maxj-2]; }
        if ( i < (maxi-1) ) { thisGamma -= vl.gammaJ()[i][maxj-2]; }
        numx += vl.endPoints()[i][maxj-1].x * fabs(thisGamma);
        numy += vl.endPoints()[i][maxj-1].y * fabs(thisGamma);
        numz += vl.endPoints()[i][maxj-1].z * fabs(thisGamma);
        den += fabs(thisGamma);
    }
    if ( den == 0.0) { 
        numx = vl.endPoints()[0][maxj-1].x; 
        numy = vl.endPoints()[0][maxj-1].y; 
        numz = vl.endPoints()[0][maxj-1].z; 
        den = 1.0;
    };
    endPoints_[0][0] = Vec3D(numx/den, numy/den, numz/den);
    gamma_[0][0] = nGamma/10.0;
    den = 0.0, numx = 0.0, numy = 0.0, numz=0.0;
    //printf("maxi of %i\n", maxi);
    //printf("maxGammaI of %i\n", maxGammaI);
    for(int i = maxGammaI; i < maxi; i++){
        thisGamma = 0.0;
        if ( i > 0 ) { thisGamma += vl.gammaJ()[i-1][maxj-2]; }
        if ( i < ( maxi-1 ) ) { thisGamma -= vl.gammaJ()[i][maxj-2]; }
        numx += vl.endPoints()[i][maxj-1].x * fabs(thisGamma);
        numy += vl.endPoints()[i][maxj-1].y * fabs(thisGamma);
        numz += vl.endPoints()[i][maxj-1].z * fabs(thisGamma);
        den += fabs(thisGamma);
    }
    if ( den == 0.0) { 
        numx = vl.endPoints()[maxi-1][maxj-1].x; 
        numy = vl.endPoints()[maxi-1][maxj-1].y; 
        numz = vl.endPoints()[maxi-1][maxj-1].z; 
        den = 1.0;
    }
    endPoints_[1][0] = Vec3D(numx/den, numy/den, numz/den);
    gamma_[1][0] = pGamma;
}

int TipFilament::ni(){
    return ni_;
}

int TipFilament::nj(){
    return nj_;
}

std::vector<std::vector<Vec3D>>& TipFilament::endPoints(){
    return endPoints_;
} 

std::vector<std::vector<Vec3D>>& TipFilament::endPointVelocity(){
    return endPointV_;
} 

std::vector<std::vector<double>>& TipFilament::gamma(){
    return gamma_;
} 

std::vector<std::vector<double>>& TipFilament::rc(){
    return rc_;
} 

std::pair<int, int> TipFilament::ijFromN( int n ){
    std::pair<int, int> ij = std::pair<int, int>();
    ij.first = n/ni_;
    ij.second = n%(ni_);

    return ij;
}

void TipFilament::initializeToHelix( Vec3D axis, double dTheta, double dZ ){
    for ( int i = 0; i < ni_; i++){
        for (int j = 1; j < nj_; j++){
            endPoints_[i][j] = endPoints_[i][j-1].rotate( Vec3D(), axis, dTheta ) - axis.norm()*dZ;
        }
    }    
}

void TipFilament::printState(){
    for (int i = 0; i < ni_; i++){
        for(int j = 0; j < nj_; j++){
            printf("loop at i = %i\n endPoint at",i);   
            endPoints_[i][j].printState();      
            endPointV_[i][j].printState();     
    
        }
    }
    for (int i = 0; i < ni_-1; i++){
        for(int j = 0; j < nj_-1; j++){
            printf("loop at i = %i rc = %5.5e \n",i, rc_[i][j]); 
        }
    }  
}
