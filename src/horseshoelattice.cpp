#include <cmath>
#include "horseshoelattice.h"
#include "vortexmath.h"
#include "stdio.h"

HorseshoeLattice::HorseshoeLattice() : ni_( 1 ), nj_( 1 ), hasTrailers_( false ), trailerVec_( 0.0, 0.0, 0.0 ), 
                                 endPoints( 2, std::vector<Vec3D>( 2, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPoints( 1, std::vector<Vec3D>( 1, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPointNormals( 1, std::vector<Vec3D>( 1, Vec3D(0.0, 0.0, 1.0) ) ), 
                                 gamma( 1, std::vector<double>( 1, 0.0  ) ) {
 
}

HorseshoeLattice::HorseshoeLattice( int ni, int nj ) : 
                                 ni_( ni ), nj_( nj ), hasTrailers_( false ), trailerVec_( 0.0, 0.0, 0.0 ),
                                 endPoints( ni+1, std::vector<Vec3D>( nj+1, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPoints( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPointNormals( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 1.0) ) ), 
                                 gamma( ni, std::vector<double>( nj, 0.0  ) ){
 
}

HorseshoeLattice::HorseshoeLattice( const HorseshoeLattice &v ) : 
                                 ni_( v.ni_ ), nj_( v.nj_ ), hasTrailers_( v.hasTrailers_ ),  trailerVec_(v.trailerVec_),
                                 endPoints( v.endPoints ),
                                 controlPoints( v.controlPoints ),
                                 controlPointNormals( v.controlPointNormals ),
                                 gamma( v.gamma ){
 
}

HorseshoeLattice::~HorseshoeLattice() {
    
}

Vec3D HorseshoeLattice::calcInfluenceCoefficient( Vec3D p, int n ){
    std::pair<int, int> ij = ijFromN( n );
    int i = ij.first; 
    int jstar = ij.second;
    Vec3D r1, r2, aOut;
    aOut = Vec3D(0.0, 0.0, 0.0); 
    //Loop through and determine influence from side legs
    for ( int j = nj_; j > jstar; j--){
        if ( j == jstar + 1) { //We're at the panel belonging to the vortex, only use 3/4 of the edge
            r1 = endPoints[i][j]  -p;
            r2 = (3.0*endPoints[i][j-1]  + 1.0*endPoints[i][j] )/4.0 - p;
            aOut += BiotSavart(r1, r2, 0.0);
        
            r1 = (3.0*endPoints[i+1][j-1] + 1.0*endPoints[i+1][j])/4.0 - p;
            r2 = endPoints[i+1][j]   - p;
            aOut += BiotSavart(r1, r2, 0.0);

        } else {
            r1 = endPoints[i][j]  -p;
            r2 = endPoints[i][j-1]-p;
            aOut += BiotSavart(r1, r2, 0.0);
        
            r1 = endPoints[i+1][j-1] - p;
            r2 = endPoints[i+1][j]   - p;
            aOut += BiotSavart(r1, r2, 0.0);
        }
    }
    if (hasTrailers_){
        r1 = endPoints[i][nj_] + trailerVec_ -p;
        r2 = endPoints[i][nj_]               -p;
        aOut += BiotSavart(r1, r2, 0.0);
        r1 = endPoints[i+1][nj_]               -p;
        r2 = endPoints[i+1][nj_] + trailerVec_ -p;
        aOut += BiotSavart(r1, r2, 0.0); 
    }
    
    //Influence of single spanwise filament
    r1 = (3.0*endPoints[i][jstar]   + 1.0*endPoints[i][jstar+1])/4.0   - p;
    r2 = (3.0*endPoints[i+1][jstar] + 1.0*endPoints[i+1][jstar+1])/4.0   - p;
    aOut += BiotSavart(r1, r2, 0.0);
    
    return aOut;
}

Vec3D HorseshoeLattice::calcInducedVelocity( Vec3D p){
    Vec3D r1, r2, aOut;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double gammaSum, gammaLocal;
    aOut = Vec3D(0.0, 0.0, 0.0);
//    #pragma omp parallel for private(gammaSum, gammaLocal, r1, r2) reduction(+:vx,vy,vz)
    for ( int i = 0; i < ni_+1; i++){
        Vec3D aTmp = Vec3D(0.0, 0.0, 0.0);//Hackery because I can't use openMP reduce on a class
        gammaSum = 0.0;
        for (int j = 0; j < nj_; j++){
            //printf("i = %i j = %i \n",i,j);
            gammaLocal = 0.0;
            if ( i < ni_ ){
                gammaLocal += gamma[i][j];
                //horseshoe segmant across quarter panel
                r1 = (3.0*endPoints[i  ][j] + 1.0*endPoints[i  ][j+1])/4.0  - p;
                r2 = (3.0*endPoints[i+1][j] + 1.0*endPoints[i+1][j+1])/4.0  - p;
                aTmp += BiotSavart(r1, r2, 0.0)*gamma[i][j];
            }
            if (i > 0) {
                gammaLocal -= gamma[i-1][j];
            }
            // Contribution of the partial horseshoe down panel edges segment
            r1 = endPoints[i][j+1] - p;
            r2 = (3.0*endPoints[i][j] + 1.0*endPoints[i][j+1])/4.0  - p;
            aTmp += BiotSavart(r1, r2, 0.0)*gammaLocal;
            
            //Contribution of complete edge filaments
            r1 = (endPoints[i][j+1] - p);
            r2 = (endPoints[i][j]   - p);
            aTmp += BiotSavart(r1, r2, 0.0)*gammaSum;
            if ( i < ni_ ){
                gammaSum   += gamma[i][j];
            }
            if (i > 0) {
                gammaSum   -= gamma[i-1][j];
            }
        }
        if (hasTrailers_){
            r1 = endPoints[i][nj_] + trailerVec_ -p;
            r2 = endPoints[i][nj_]               -p;
            aTmp += BiotSavart(r1, r2, 0.0) * gammaSum;
        }
        vx += aTmp.x;
        vy += aTmp.y;
        vz += aTmp.z;
        
    }
    aOut.x = vx; aOut.y = vy; aOut.z = vz;
    return aOut;
}

void HorseshoeLattice::snapToUnit(){
    Vec3D dx = Vec3D(1.0/(double)( ni_ ), 0.0                , 0.0);
    Vec3D dy = Vec3D(0.0                , 1.0/(double)( nj_ ), 0.0);
    
    for (int i = 0; i < ni_+1; i++){
        if (i > 0){
            endPoints[i][0] = endPoints[i-1][0] + dx;
        } else {
            endPoints[i][0] = Vec3D(0.0, 0.0, 0.0);
        }   
        for (int j = 1; j < nj_+1; j++){
           endPoints[i][j] = endPoints[i][j-1] + dy;
        }
    }
    centerControlPoints();
}

void HorseshoeLattice::snapToAspectTaper( double ar, double taper ){
    snapToAspectTaperSweep( ar, taper, 0.0 );
}

void HorseshoeLattice::snapToAspectTaperSweep( double ar, double taper, double sweep ){
    double cr = 1.0, ct = cr*taper, cbar = (cr+ct)/2.0, b = ar*cbar, cLocal; //s = b*cbar
    Vec3D dx = Vec3D(b/(double)( ni_ ), 0.0                , 0.0);
    Vec3D dy = Vec3D(0.0                , cr/(double)( nj_ ), 0.0);
    double yLe;
    for (int i = 0; i < ni_+1; i++){
        if (i > 0){
            endPoints[i][0] = endPoints[i-1][0] + dx;
        } else {
            endPoints[i][0] = Vec3D(0.0, 0.0, 0.0);
        }
        cLocal = cr - (cr-ct)*endPoints[i][0].x/b;
        yLe = -cLocal/4.0 + endPoints[i][0].x*tan(sweep);
        dy.y = cLocal/(double)( nj_ );
        endPoints[i][0].y = yLe;
        for (int j = 1; j < nj_+1; j++){
           endPoints[i][j] = endPoints[i][j-1] + dy;
        }
    }
    centerControlPoints();
}

void HorseshoeLattice::flipTip( double dihedralBreak, double dihedral ) {
    Vec3D leVector = (dihedralBreak*endPoints[ni_][0  ] + (1.0-dihedralBreak)*endPoints[0][0  ] ); 
    Vec3D teVector = (dihedralBreak*endPoints[ni_][nj_] + (1.0-dihedralBreak)*endPoints[0][nj_] ); 
    Vec3D rotateAxis = leVector-teVector;
    
    for (int i = 0; i < ni_+1; i++){
        Vec3D localVec = endPoints[i][0] - endPoints[0][0]; 
        if ( localVec.magnitude() < leVector.magnitude() ) continue;
        for (int j = 0; j < nj_+1; j++){
            endPoints[i][j] = endPoints[i][j].rotate(leVector, rotateAxis, dihedral);
        }
    }   
            
    centerControlPoints();
}

void HorseshoeLattice::centerControlPoints(){
    for (int i = 0; i < ni_; i++){
        for (int j = 0; j < nj_; j++){
            controlPoints[i][j] = ( 1.0*endPoints[i  ][j  ] +
                                    1.0*endPoints[i+1][j  ] +
                                    3.0*endPoints[i  ][j+1] +
                                    3.0*endPoints[i+1][j+1] )/ 8.0;
        }
    }
    calcControlPointNormals();
}

void HorseshoeLattice::calcControlPointNormals(){
    Vec3D v1, v2, n1, n2;
    for (int i = 0; i < ni_; i++){
        for (int j = 0; j < nj_; j++){
            v1 = endPoints[i  ][j+1] - endPoints[i  ][j  ];
            v2 = endPoints[i+1][j  ] - endPoints[i  ][j  ];
            n1 = v1.cross(v2).norm();
            v1 = endPoints[i+1][j+1] - endPoints[i  ][j+1];
            v2 = endPoints[i+1][j+1] - endPoints[i+1][j  ];
            n2 = v2.cross(v1).norm();
            controlPointNormals[i][j] = (n1+n2)/2.0;
        }
    } 
}

double HorseshoeLattice::dSpan(int i, int j){
    //dSpan is length of vector
    Vec3D v1 = endPoints[i  ][j] - endPoints[i][j+1];
    Vec3D v2 = endPoints[i+1][j] - endPoints[i][j+1];
    Vec3D vn = v2-(v1.norm()*v2.dot(v1.norm()));
    return  vn.magnitude();
}

Vec3D HorseshoeLattice::gammaVector( int i, int j ){
    Vec3D vOut = gamma[i][j] * (   
        (3.0*endPoints[i+1][j  ] +
         1.0*endPoints[i+1][j+1])  -
        (3.0*endPoints[i  ][j  ] +
         1.0*endPoints[i  ][j+1]) ) / 4.0;
    return vOut; 
}

Vec3D HorseshoeLattice::gammaCenterPoint( int i, int j ){
    Vec3D vOut = (   
        (3.0*endPoints[i+1][j  ] +
         1.0*endPoints[i+1][j+1])  +
        (3.0*endPoints[i  ][j  ] +
         1.0*endPoints[i  ][j+1]) ) / 8.0;
    return vOut; 
}

void HorseshoeLattice::rotate( Vec3D point, Vec3D axis, double theta ){
    for (int i = 0; i < ni_; i++){
        for (int j = 0; j < nj_; j++){
            controlPoints[i][j] = controlPoints[i][j].rotate(point, axis, theta);
            controlPointNormals[i][j] = controlPointNormals[i][j].rotate(Vec3D(0.0, 0.0, 0.0), axis, theta);
        }
    }
     
    for (int i = 0; i < ni_+1; i++){
        for (int j = 0; j < nj_+1; j++){
            endPoints[i][j] = endPoints[i][j].rotate(point, axis, theta);
        }
    }
}

void HorseshoeLattice::translate( Vec3D dir ){
    for (int i = 0; i < ni_; i++){
        for (int j = 0; j < nj_; j++){
            controlPoints[i][j] += dir;
        }
    }
    for (int i = 0; i < ni_+1; i++){
        for (int j = 0; j < nj_+1; j++){
            endPoints[i][j] += dir;
        }
    }
}   

int HorseshoeLattice::ni(){
    return ni_;
}

int HorseshoeLattice::nj(){
    return nj_;
}

std::vector<std::vector<Vec3D>>& HorseshoeLattice::getEndPoints(){
    return endPoints;
} 

std::vector<std::vector<Vec3D>>& HorseshoeLattice::getControlPoints(){
    return controlPoints;
} 

std::vector<std::vector<Vec3D>>& HorseshoeLattice::getControlPointNormals(){
    return controlPointNormals;
}
 
std::vector<std::vector<double>>& HorseshoeLattice::getGamma(){
    return gamma;
} 

void HorseshoeLattice::setHasTrailers( bool b ) {
    hasTrailers_ = b;
}

bool HorseshoeLattice::hasTrailers(){
    return hasTrailers_;
}

Vec3D HorseshoeLattice::trailerVec(){
    return trailerVec_;
}

void HorseshoeLattice::setTrailerVec( Vec3D v ){
   trailerVec_ = v;
}
    
   
std::pair<int, int> HorseshoeLattice::ijFromN( int n ){
    std::pair<int, int> ij = std::pair<int, int>();
    ij.first  = n/(nj_);
    ij.second = n%(nj_);

    return ij;
}

void HorseshoeLattice::printState(){
    for (int i = 0; i < ni_+1; i++){
        for(int j = 0; j < nj_+1; j++){
            printf("loop at i = %i\n endPoint at",i);   
            endPoints[i][j].printState();      
        }
    }
    for (int i = 0; i < ni_; i++){
        for(int j = 0; j < nj_; j++){
            printf("loop at i = %i\n controlPoint at",i);   
            controlPoints[i][j].printState();      
        }
    }
}

int HorseshoeLattice::ijToN(int i, int j){
    int n = i*(nj_)+j;
    return n;
}

int HorseshoeLattice::maxN(){
    return ni_*nj_;
}
