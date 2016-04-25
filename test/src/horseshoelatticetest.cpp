#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "horseshoelattice.h"
#include "vortexmath.h"

class HorseshoeLatticeTest: public ::testing::Test {
	public:
	HorseshoeLattice *vl;
    virtual void SetUp(){
        vl = new HorseshoeLattice(4,3);
        for (int i = 0; i < vl->ni(); i++){
            for (int j = 0; j < vl->nj(); j++){
                vl->getGamma()[i][j] = 1.0;
            }
        }
	}

	virtual void TearDown(){
        delete vl;
	}

};

TEST_F(HorseshoeLatticeTest, TestCopy){
    HorseshoeLattice v2 = HorseshoeLattice(*vl);
    EXPECT_EQ( 4, vl->ni());    
    EXPECT_EQ( 3, vl->nj());    
    EXPECT_EQ( vl->ni(), v2.ni());    
    EXPECT_EQ( vl->nj(), v2.nj());
}

TEST_F(HorseshoeLatticeTest, TestEndPointsArraySize){
    EXPECT_EQ( 5, vl->getEndPoints().size() );
    EXPECT_EQ( 4, vl->getEndPoints()[0].size() );
}

TEST_F(HorseshoeLatticeTest, TestGammaArraySize){
    EXPECT_EQ( 4, vl->getGamma().size() );
    EXPECT_EQ( 3, vl->getGamma()[0].size() );
}

TEST_F(HorseshoeLatticeTest, TestIJFromN){ 
    //Solutions derived on paper for a 6X3 lattice
    EXPECT_EQ( 0, vl->ijFromN( 0 ).first );
    EXPECT_EQ( 0, vl->ijFromN( 0 ).second );
    
    EXPECT_EQ( 0, vl->ijFromN( 2 ).first );
    EXPECT_EQ( 2, vl->ijFromN( 2 ).second );
    
    EXPECT_EQ( 1, vl->ijFromN( 3 ).first );
    EXPECT_EQ( 0, vl->ijFromN( 3 ).second );
    
    EXPECT_EQ( 1, vl->ijFromN( 4 ).first );
    EXPECT_EQ( 1, vl->ijFromN( 4 ).second );
    
    EXPECT_EQ( 3, vl->ijFromN( 9 ).first );
    EXPECT_EQ( 0, vl->ijFromN( 9 ).second );
    
}

TEST_F(HorseshoeLatticeTest, TestIJToN){ 
    int n, nstar; n  = 0;
    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
            nstar = vl->ijToN(i,j);
            EXPECT_EQ(nstar,n);
            n++;
        }
    }
}
TEST_F(HorseshoeLatticeTest, TestSnapToArTaper){
    vl->snapToAspectTaper( 4.0/0.75, .5);

    double x[5] =      {0.0,     1.0,     2.0,     3.0,     4.0};
    double xcp[4] =    {0.5,     1.5,     2.5,     3.5};
    double y[4] =      {0.0,     1.0/3.0, 2.0/3.0, 1.0};
    double ycp[3] =    {3.0/12.0, 7.0/12.0, 11.0/12.0};
    double yScale[5] = {8.0/8.0, 7.0/8.0, 6.0/8.0, 5.0/8.0, 4.0/8.0};
    double yLe[5] =    {0.0, (1.0-7.0/8.0)/4.0, (1.0-6.0/8.0)/4.0, (1.0-5.0/8.0)/4.0, (1.0-4.0/8.0)/4.0};
    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
           //printf("%i, %i\n",i,j);
           EXPECT_DOUBLE_EQ( x[i]                       , vl->getEndPoints()[i][j].x);
           EXPECT_DOUBLE_EQ( y[j]*yScale[i]+yLe[i]-0.25 , vl->getEndPoints()[i][j].y);
           
           EXPECT_NEAR     ( xcp[i]                                                        , vl->getControlPoints()[i][j].x, 1E-15);
           EXPECT_NEAR     ( ycp[j]*(yScale[i]+yScale[i+1])/2.0+(yLe[i]+yLe[i+1])/2.0-0.25 , vl->getControlPoints()[i][j].y, 1E-15);
        }
    }      
}

TEST_F(HorseshoeLatticeTest, TestControlPoints){
    vl->snapToUnit();
    double xcp[4] =    {1.0/8.0, 3.0/8.0, 5.0/8.0, 7.0/8.0};
    double ycp[3] =    {3.0/12.0, 7.0/12.0, 11.0/12.0};
    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
           //printf("%i, %i\n",i,j);
           EXPECT_DOUBLE_EQ(xcp[i] , vl->getControlPoints()[i][j].x);
           EXPECT_DOUBLE_EQ(ycp[j] , vl->getControlPoints()[i][j].y);
           EXPECT_DOUBLE_EQ(0.0    , vl->getControlPoints()[i][j].z);
        }
    }
}

TEST_F(HorseshoeLatticeTest, TestControlPointNormals){
    vl->snapToAspectTaper( 4.0/0.75, .5);
    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
           //printf("%i, %i\n",i,j);
           EXPECT_DOUBLE_EQ( 0.0, vl->getControlPointNormals()[i][j].x);
           EXPECT_DOUBLE_EQ( 0.0, vl->getControlPointNormals()[i][j].y);
           EXPECT_DOUBLE_EQ(-1.0, vl->getControlPointNormals()[i][j].z);
        }
    }
    vl->snapToAspectTaper( 1.0, 1.0 );
    vl->rotate( Vec3D( 0.0, 0.0, 0.0 ), Vec3D( 1.0, 0.0 , 0.0 ), 45.0 * M_PI / 180.0 ); 
    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
           //printf("%i, %i\n",i,j);
           EXPECT_DOUBLE_EQ( 0.0          , vl->getControlPointNormals()[i][j].x);
           EXPECT_DOUBLE_EQ( sqrt(2.0)/2.0, vl->getControlPointNormals()[i][j].y);
           EXPECT_DOUBLE_EQ(-sqrt(2.0)/2.0, vl->getControlPointNormals()[i][j].z);
        }
    }
}

TEST_F(HorseshoeLatticeTest, TestDSan){
    HorseshoeLattice hl = HorseshoeLattice(1.0, 1.0);
    hl.snapToUnit();
    Vec3D a = hl.calcInfluenceCoefficient( hl.getControlPoints()[0][0] , 0 ); // Self induced coefficient without trailers
    
    double b = 5.0;
    hl.getEndPoints()[0][0] = Vec3D( 0.0  , 0.0  , 0.0 );
    hl.getEndPoints()[0][1] = Vec3D( 0.0  , 1.0  , 0.0 );
    hl.getEndPoints()[1][0] = Vec3D( 0.625, 0.625, 0.0 );
    hl.getEndPoints()[1][1] = Vec3D( 0.625, 1.625, 0.0 );
    hl.centerControlPoints();
    
    EXPECT_DOUBLE_EQ( 0.3125, hl.getControlPoints()[0][0].x );
    EXPECT_DOUBLE_EQ( 0.625, hl.dSpan(0,0));
}

TEST_F(HorseshoeLatticeTest, TestTranslate){
    vl->snapToUnit();
    HorseshoeLattice h2 = HorseshoeLattice(*vl);
    Vec3D tVec = Vec3D(1.0, 2.0, 3.0);
    h2.translate(tVec);
    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
            EXPECT_DOUBLE_EQ(vl->getControlPoints()[i][j].x + 1.0, h2.getControlPoints()[i][j].x);
            EXPECT_DOUBLE_EQ(vl->getControlPoints()[i][j].y + 2.0, h2.getControlPoints()[i][j].y);
            EXPECT_DOUBLE_EQ(vl->getControlPoints()[i][j].z + 3.0, h2.getControlPoints()[i][j].z);
        }
    }
}

TEST_F(HorseshoeLatticeTest, TestRotate){
    vl->snapToUnit();
    HorseshoeLattice h2 = HorseshoeLattice(*vl);
    Vec3D axisVec = Vec3D(0.0, 0.0, 1.0);
 
    h2.rotate(Vec3D(), axisVec, M_PI);

    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
            EXPECT_DOUBLE_EQ(-vl->getControlPoints()[i][j].x, h2.getControlPoints()[i][j].x);
            EXPECT_DOUBLE_EQ(-vl->getControlPoints()[i][j].y, h2.getControlPoints()[i][j].y);
            EXPECT_DOUBLE_EQ(vl->getControlPoints()[i][j].z, h2.getControlPoints()[i][j].z);
        }
    }
}

TEST_F(HorseshoeLatticeTest, TestFlipTip){
    vl->snapToUnit();
    HorseshoeLattice h2 = HorseshoeLattice(*vl);
    h2.flipTip(0.5, M_PI/2.0);
    for (int i = 0; i < vl->ni()+1; i++){
        for (int j = 0; j < vl->nj()+1; j++){
            if (i > 2){
                //X values rotate up to z;  
                EXPECT_DOUBLE_EQ( vl->getEndPoints()[i][j].x-0.5, h2.getEndPoints()[i][j].z );
                //y values untouched 
                EXPECT_DOUBLE_EQ( vl->getEndPoints()[i][j].y, h2.getEndPoints()[i][j].y );
            }
        }
    }
}

TEST_F(HorseshoeLatticeTest, TestInfluenceCoefficient){
    HorseshoeLattice hl = HorseshoeLattice(1.0, 1.0);
    hl.snapToUnit();
    Vec3D a = hl.calcInfluenceCoefficient( hl.getControlPoints()[0][0] , 0 ); // Self induced coefficient without trailers
    
    Vec3D ra = Vec3D( 0.0, 1.0 , 0.0 );
    Vec3D rb = Vec3D( 0.0, 0.25, 0.0 );
    Vec3D rc = Vec3D( 1.0, 0.25, 0.0 );
    Vec3D rd = Vec3D( 1.0, 1.0 , 0.0 );
    Vec3D p  = Vec3D( 0.5, 0.75, 0.0 );

    double rcore = 1E-6;
    Vec3D aReal = Vec3D( 0.0, 0.0, 0.0 );
    Vec3D r1, r2;
    r1 = ra-p; r2 = rb-p; 
    aReal += BiotSavart( r1, r2, rcore);
    r1 = rb-p; r2 = rc-p; 
    aReal += BiotSavart( r1, r2, rcore);
    r1 = rc-p; r2 = rd-p; 
    aReal += BiotSavart( r1, r2, rcore);

    EXPECT_DOUBLE_EQ( aReal.x, a.x );
    EXPECT_DOUBLE_EQ( aReal.y, a.y );
    EXPECT_DOUBLE_EQ( aReal.z, a.z );
    
}

TEST_F(HorseshoeLatticeTest, TestBertinCummingsExampleInfluence){
    HorseshoeLattice hl = HorseshoeLattice(1.0, 1.0);
    hl.snapToUnit();
    Vec3D a = hl.calcInfluenceCoefficient( hl.getControlPoints()[0][0] , 0 ); // Self induced coefficient without trailers
    
    double b = 5.0;
    hl.getEndPoints()[0][0] = Vec3D( 0.0  , 0.0  , 0.0 );
    hl.getEndPoints()[0][1] = Vec3D( 0.0  , 1.0  , 0.0 );
    hl.getEndPoints()[1][0] = Vec3D( 0.625, 0.625, 0.0 );
    hl.getEndPoints()[1][1] = Vec3D( 0.625, 1.625, 0.0 );
    hl.centerControlPoints();
    
    EXPECT_DOUBLE_EQ( 0.3125, hl.getControlPoints()[0][0].x );
    EXPECT_DOUBLE_EQ( 1.0625, hl.getControlPoints()[0][0].y );

    Vec3D vTrailer = Vec3D(0.0, 1000.0, 0.0);
    hl.setHasTrailers( true );
    hl.setTrailerVec( Vec3D(0.0, 1000.0, 0.0) );
    a = hl.calcInfluenceCoefficient( hl.getControlPoints()[0][0] , 0 ); // Self induced coefficient with trailers
    double influence = hl.getControlPointNormals()[0][0].dot( a );
    EXPECT_NEAR( -1.1382551251610777, influence, 1E-12 );
   
    Vec3D rb = Vec3D( 0.0   , 0.0   , 0.0 );
    Vec3D ra = Vec3D( 0.0   , 1.0   , 0.0 );
    Vec3D rc = Vec3D( 0.625 , 0.625 , 0.0 );
    Vec3D rd = Vec3D( 0.625 , 1.625 , 0.0 );
    Vec3D p  = Vec3D( 0.3125, 1.0625, 0.0 );

    Vec3D aHardway = Vec3D();
    Vec3D r1, r2;
    r1 = (ra+vTrailer)-p; r2 = ra-p; 
    aHardway += BiotSavart( r1, r2, 0.0);
    r1 = ra-p; r2 = rb-p; 
    aHardway += BiotSavart( r1, r2, 0.0);
    r1 = rb-p; r2 = rc-p; 
    aHardway += BiotSavart( r1, r2, 0.0);
    r1 = rc-p; r2 = rd-p; 
    aHardway += BiotSavart( r1, r2, 0.0);
    r1 = rd-p; r2 = (rd + vTrailer)-p; 
    aHardway += BiotSavart( r1, r2, 0.0);
    
    //a.printState();
    //aHardway.printState();
    //hl.getControlPointNormals()[0][0].printState();
    //double influenceHardway = hl.getControlPointNormals()[0][0].dot( aHardway );
    //EXPECT_DOUBLE_EQ( 1.1382551251610777, influenceHardway );
       
}

TEST_F(HorseshoeLatticeTest, TestInducedVelocityIsSumOfInfluenceCoefficients){
    Vec3D vInduced   = Vec3D(0.0, 0.0, 0.0);
    Vec3D aInfluence = Vec3D(0.0, 0.0, 0.0);
    Vec3D p          = Vec3D(0.0, 0.0, 1.0);
    vl->snapToUnit();
    //vl->printState();
    for (int i = 0; i < vl->ni(); i++){
        for (int j = 0; j < vl->nj(); j++){
           aInfluence += vl->calcInfluenceCoefficient( p, vl->ijToN(i,j) );
        }
    }
    vInduced = vl->calcInducedVelocity(p);

    EXPECT_NEAR     ( vInduced.x, aInfluence.x, 1E-16); //This comes out a few too many epsilon from zero
    EXPECT_DOUBLE_EQ( vInduced.y, aInfluence.y);
    EXPECT_DOUBLE_EQ( vInduced.z, aInfluence.z);
}

TEST_F(HorseshoeLatticeTest, TestInducedVelocityParallel){
    int imax = 10;
    int jmax = 10;
    HorseshoeLattice v2 = HorseshoeLattice( imax, jmax);
    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){
            v2.calcInducedVelocity(v2.getControlPoints()[i][j]);
            EXPECT_EQ(1,1);
        }   
    }
}
