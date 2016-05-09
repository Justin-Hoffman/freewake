#include <cmath>
#include "horseshoelattice.h"
#include "vortexmath.h"
#include "stdio.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

HorseshoeLattice::HorseshoeLattice() : rc_(1E-6), ni_( 1 ), nj_( 1 ), chordwiseSpacing_( PointSpacing::Linear ), spanwiseSpacing_(PointSpacing::Linear), 
                                 hasTrailers_( false ), trailerVec_( 0.0, 0.0, 0.0 ), 
                                 endPoints( 2, std::vector<Vec3D>( 2, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPoints( 1, std::vector<Vec3D>( 1, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPointNormals( 1, std::vector<Vec3D>( 1, Vec3D(0.0, 0.0, 1.0) ) ), 
                                 gamma( 1, std::vector<double>( 1, 0.0  ) ), 
                                 trailedPoints_( ni_+1, Vec3D()  ),
                                 netGamma_( 1, 0.0 ){
 
}

HorseshoeLattice::HorseshoeLattice( int ni, int nj ) : 
                                 rc_(1E-6), ni_( ni ), nj_( nj ), chordwiseSpacing_( PointSpacing::Linear ), spanwiseSpacing_(PointSpacing::Linear), 
                                 hasTrailers_( false ), trailerVec_( 0.0, 0.0, 0.0 ),
                                 endPoints( ni+1, std::vector<Vec3D>( nj+1, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPoints( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 0.0) ) ), 
                                 controlPointNormals( ni, std::vector<Vec3D>( nj, Vec3D(0.0, 0.0, 1.0) ) ), 
                                 gamma( ni, std::vector<double>( nj, 0.0  ) ),
                                 trailedPoints_( ni_+1, Vec3D() ),
                                 netGamma_(ni, 0.0){
 
}

HorseshoeLattice::HorseshoeLattice( const HorseshoeLattice &v ) : 
                                 rc_(v.rc_),ni_( v.ni_ ), nj_( v.nj_ ), chordwiseSpacing_( v.chordwiseSpacing_ ), spanwiseSpacing_( v.spanwiseSpacing_ ), 
                                 hasTrailers_( v.hasTrailers_ ),  trailerVec_(v.trailerVec_),
                                 endPoints( v.endPoints ),
                                 controlPoints( v.controlPoints ),
                                 controlPointNormals( v.controlPointNormals ),
                                 gamma( v.gamma ),
                                 trailedPoints_( v.trailedPoints_ ),
                                 netGamma_( v.netGamma_ ){
 
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
            aOut += BiotSavart(r1, r2, rc_);
        
            r1 = (3.0*endPoints[i+1][j-1] + 1.0*endPoints[i+1][j])/4.0 - p;
            r2 = endPoints[i+1][j]   - p;
            aOut += BiotSavart(r1, r2, rc_);

        } else {
            r1 = endPoints[i][j]  -p;
            r2 = endPoints[i][j-1]-p;
            aOut += BiotSavart(r1, r2, rc_);
        
            r1 = endPoints[i+1][j-1] - p;
            r2 = endPoints[i+1][j]   - p;
            aOut += BiotSavart(r1, r2, rc_);
        }
    }
    if (hasTrailers_){
        r1 = endPoints[i][nj_] + trailerVec_ -p;
        r2 = endPoints[i][nj_]               -p;
        aOut += BiotSavart(r1, r2, rc_);
        r1 = endPoints[i+1][nj_]               -p;
        r2 = endPoints[i+1][nj_] + trailerVec_ -p;
        aOut += BiotSavart(r1, r2, rc_); 
    } else {
        r1 = trailedPoints_[i]-p;
        r2 = endPoints[i][nj_]-p;
        aOut += BiotSavart(r1, r2, rc_);
        r1 = endPoints[i+1][nj_]-p;
        r2 = trailedPoints_[i+1]-p;
        aOut += BiotSavart(r1, r2, rc_); 
    } 
    
    //Influence of single spanwise filament
    r1 = (3.0*endPoints[i][jstar]   + 1.0*endPoints[i][jstar+1])/4.0   - p;
    r2 = (3.0*endPoints[i+1][jstar] + 1.0*endPoints[i+1][jstar+1])/4.0   - p;
    aOut += BiotSavart(r1, r2, rc_);
    
    return aOut;
}

Vec3D HorseshoeLattice::calcInducedVelocity( Vec3D p, int jStart ){
    Vec3D r1, r2, aOut;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double gammaSum, gammaLocal;
    aOut = Vec3D(0.0, 0.0, 0.0);
    #pragma omp parallel for private(gammaSum, gammaLocal, r1, r2) reduction(+:vx,vy,vz)
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
                aTmp += BiotSavart(r1, r2, rc_)*gamma[i][j];
            }
            if (i > 0) {
                gammaLocal -= gamma[i-1][j];
            }
            // Contribution of the partial horseshoe down panel edges segment
            r1 = endPoints[i][j+1] - p;
            r2 = (3.0*endPoints[i][j] + 1.0*endPoints[i][j+1])/4.0  - p;
            aTmp += BiotSavart(r1, r2, rc_)*gammaLocal;
            
            //Contribution of complete edge filaments
            r1 = (endPoints[i][j+1] - p);
            r2 = (endPoints[i][j]   - p);
            aTmp += BiotSavart(r1, r2, rc_)*gammaSum;
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
            aTmp += BiotSavart(r1, r2, rc_) * gammaSum;
        } else {
            double dx = 1.0E12;;
            if (i < ni_) dx = fmin(dx,(endPoints[i][nj_]-endPoints[i+1][nj_]).magnitude());
            if (i > 0  ) dx = fmin(dx,(endPoints[i][nj_]-endPoints[i-1][nj_]).magnitude());
            r1 = trailedPoints_[i]-p;
            r2 = endPoints[i][nj_]-p;
            aTmp += BiotSavart(r1, r2, fmax(dx,fmax(rc_,dx)) ) * gammaSum;
    
        }
        vx += aTmp.x;
        vy += aTmp.y;
        vz += aTmp.z;
        
    }
    aOut.x = vx; aOut.y = vy; aOut.z = vz;
    return aOut;
}

void HorseshoeLattice::snapToUnit(){
    Vec3D dy = Vec3D(0.0                , 1.0/(double)( ni_ ), 0.0);
    Vec3D dx = Vec3D(-1.0/(double)( nj_ ), 0.0                , 0.0);
    
    for (int i = 0; i < ni_+1; i++){
        if (i > 0){
            endPoints[i][0] = endPoints[i-1][0] + dy;
        } else {
            endPoints[i][0] = Vec3D(0.0, 0.0, 0.0);
        }   
        for (int j = 1; j < nj_+1; j++){
           endPoints[i][j] = endPoints[i][j-1] + dx;
        }
        trailedPoints_[i] = endPoints[i][nj_];
    }
    centerControlPoints();
}

void HorseshoeLattice::snapToAspectTaper( double ar, double taper ){
    snapToAspectTaperSweep( ar, taper, 0.0 );
}

void HorseshoeLattice::snapToAspectTaperSweep( double ar, double taper, double sweep ){
    double cr = 1.0, ct = cr*taper, cbar = (cr+ct)/2.0, b = ar*cbar, cLocal; //s = b*cbar
    
    std::vector<double> span_loc = std::vector<double>(ni_+1, 0.0);
    std::vector<double> chord_loc = std::vector<double>(nj_+1, 0.0);
    switch (spanwiseSpacing_){
        case PointSpacing::Linear :
            for(int i = 1; i < ni_+1; i++){ 
                span_loc[i] = (double) i * b / ((double) ni_); 
            }
            break;
        case PointSpacing::Cosine :
            for(int i = 1; i < ni_+1; i++){ 
                span_loc[i] = ( cos( M_PI -  (double) i * M_PI / ((double) ni_) ) + 1.0 ) * b/2.0;
            }
            break;
        case PointSpacing::HalfCosine :
            for(int i = 1; i < ni_+1; i++){ 
                span_loc[i] = ( cos( M_PI/2.0 - (double) i * M_PI / (2.0 * (double) ni_) ) ) * b;
            }
            break;
    }
    switch( chordwiseSpacing_ ){
        case PointSpacing::Linear :
            for(int j = 1; j < nj_+1; j++){
                chord_loc[j] = (double) j * 1.0 / ((double) nj_);
            }
            break;
        case PointSpacing::Cosine :
            for(int j = 1; j < nj_+1; j++){ 
                chord_loc[j] = ( cos( M_PI -  (double) j * M_PI / ((double) nj_) ) + 1.0 ) * 1.0/2.0;
            }
            break;
        case PointSpacing::HalfCosine :
            for(int j = 1; j < nj_+1; j++){ 
                chord_loc[j] = ( cos( M_PI/2.0 - (double) j * M_PI / (2.0 * (double) nj_) ) );
            }
            break;
    }


    Vec3D dy = Vec3D( 0.0                , b/(double)( ni_ ), 0.0);
    Vec3D dx = Vec3D(-1.0/(double)( nj_ ), 0.0              , 0.0);
    double xLe;
    for (int i = 0; i < ni_+1; i++){
        endPoints[i][0] = Vec3D(0.0, span_loc[i], 0.0);
        cLocal = cr - (cr-ct)*span_loc[i]/b;
        xLe =  cLocal/4.0 - span_loc[i]*tan(sweep);
        endPoints[i][0].x = xLe;
        for (int j = 1; j < nj_+1; j++){
           endPoints[i][j].x = xLe-chord_loc[j]*cLocal;
           endPoints[i][j].y = span_loc[i];
           endPoints[i][j].z =  0.0;
        }
        trailedPoints_[i] = endPoints[i][nj_];
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

double HorseshoeLattice::dChord(int i ){
    //dSpan is length of vector
    Vec3D v1 = endPoints[i][nj_-1] - endPoints[i][0];
    return  v1.magnitude();
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

void HorseshoeLattice::scale( double sc ){
    for (int i = 0; i < ni_; i++){
        for (int j = 0; j < nj_; j++){
            controlPoints[i][j] *= sc;
        }
    }
    for (int i = 0; i < ni_+1; i++){
        for (int j = 0; j < nj_+1; j++){
            endPoints[i][j] *= sc;
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

double HorseshoeLattice::getRc(){
    return rc_;
}

void HorseshoeLattice::setRc( double rc){
    rc_ = rc;
}


PointSpacing HorseshoeLattice::chordwiseSpacing(){
    return chordwiseSpacing_;
}

PointSpacing HorseshoeLattice::spanwiseSpacing(){
    return spanwiseSpacing_;
}

void HorseshoeLattice::chordwiseSpacing( PointSpacing ps ){
    chordwiseSpacing_ = ps;
}

void HorseshoeLattice::spanwiseSpacing( PointSpacing ps ){
    spanwiseSpacing_ = ps;
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

std::vector<Vec3D>& HorseshoeLattice::getTrailedPoints(){
    return trailedPoints_;
} 

std::vector<double>& HorseshoeLattice::getNetGamma(){
    return netGamma_;
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
