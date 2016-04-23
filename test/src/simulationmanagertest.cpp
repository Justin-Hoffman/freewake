#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
    
#include "simulationmanager.h"

class SimulationManagerTest: public ::testing::Test {
	public:
    SimulationManager* sm;
	virtual void SetUp(){
       sm = new SimulationManager();
	}

	virtual void TearDown(){
        delete sm;
	}

};

TEST_F(SimulationManagerTest, TestLiftingLine){
    LiftingSurface ls = LiftingSurface(10,1);
    ls.setAspectRatio( 10.0 );
    ls.setPitch( 0.0 * M_PI / 180.0 );
    ls.getLattice().setHasTrailers(true);   
    ls.getLattice().setTrailerVec( Vec3D(0.0, 1000.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(1.0, 0.0, 0.0), -10.0*M_PI/180.0 ) );   
    ls.updateLattice( );
    ls.getLattice().printState();
    sm->addSurface(&ls);
    sm->setGlobalLinearVelocity( Vec3D(0.0, 100.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(1.0,0.0, 0.0), -10.0*M_PI/180.0 ) );
    sm->solve();   
    double gamma[10] = { 35.03083712393401328, 43.13601948056102486, 46.29353105076378938, 47.71129676895412075, 48.28231926132752250, 48.28231926132752960, 47.71129676895412075, 46.29353105076377517, 43.13601948056103197, 35.03083712393401328};

    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_DOUBLE_EQ( gamma[n], ls.getLattice().getGamma()[i][j] );
        }
    }
     
}

TEST_F(SimulationManagerTest, TestLiftingLineTwoChord){
    LiftingSurface ls = LiftingSurface(10,2);
    ls.setAspectRatio( 10.0 );
    ls.setPitch( 10.0 * M_PI / 180.0 );
    ls.getLattice().setHasTrailers(true);   
    ls.getLattice().setTrailerVec( Vec3D(0.0, 1000.0, 0.0) );   
    ls.updateLattice( );
    sm->addSurface(&ls);
    sm->setGlobalLinearVelocity( Vec3D(0.0, 100.0, 0.0) );
    sm->solve();   
    double gamma[10] = { 35.03083712393401328, 43.13601948056102486, 46.29353105076378938, 47.71129676895412075, 48.28231926132752250, 48.28231926132752960, 47.71129676895412075, 46.29353105076377517, 43.13601948056103197, 35.03083712393401328};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        double gammaSum = 0.0;
        n++;
        for(int j = 0; j < ls.nChord(); j++){
            gammaSum += ls.getLattice().getGamma()[i][j];
        }    
        EXPECT_NEAR( gamma[n], gammaSum, gamma[n]*0.01); //Two chordwise station should find net circulation to within 1% of one panel
    }
     
}

TEST_F(SimulationManagerTest, TestBertinCummingsExampleSweptWing){
    LiftingSurface ls = LiftingSurface(4,1);
    ls.setAspectRatio(2.5);
    ls.setSweep( 45.0 * M_PI / 180.0 );
    ls.getLattice().setHasTrailers(true);   
    ls.getLattice().setTrailerVec( Vec3D(0.0, 1000.0, 0.0).rotate( Vec3D(0.0,0.0,0.0), Vec3D(1.0,0.0,0.0), -1.0*M_PI/180.0) );   
    ls.updateLattice( );
    ls.getLattice().printState();
    sm->addSurface(&ls);

    LiftingSurface ls2 = LiftingSurface( ls );
    ls2.setSweep( -45.0 * M_PI / 180.0 );
    ls2.updateLattice( );
    ls2.getLattice().translate( Vec3D(-2.5, 2.5, 0.0) );
    ls2.getLattice().printState();
    sm->addSurface(&ls2);
   
    EXPECT_TRUE( ls2.getLattice().hasTrailers() ); 
    sm->setGlobalLinearVelocity( Vec3D(0.0, 100.0, 0.0).rotate( Vec3D(0,0,0), Vec3D(1,0,0), -1.0*M_PI/180.0) );
    sm->solve();   
    //Taken from 7.49a/b/c/d in sixth edition of Bertin and Cummings - Note the data is given only to 3 sig figs
    double gammaBertinCummings[4] = {2.993628011055486, 3.147147396237819, 3.136181725867653, 2.741417592541655};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_NEAR( gammaBertinCummings[n], ls.getLattice().getGamma()[i][j], 5E-3 );
        }
    }      
    for(int i = 0; i < ls2.nSpan(); i++){
        for(int j = 0; j < ls2.nChord(); j++){
            EXPECT_NEAR( gammaBertinCummings[n], ls2.getLattice().getGamma()[i][j], 5E-3 );
            n--;
        }
    }
    EXPECT_NEAR( 0.0601, sm->netLift()/(1.0/2.0 * 1.0 * 100.0 * 100.0 * 5.0 ), 5E-5);
   
    sm->setGlobalLinearVelocity( Vec3D(0.0, 100.0, 0.0).rotate( Vec3D(0,0,0), Vec3D(1,0,0), -2.0*M_PI/180.0) );
    sm->solve();   
    EXPECT_NEAR( 0.1202, sm->netLift()/(1.0/2.0 * 1.0 * 100.0 * 100.0 * 5.0 ), 1E-4);
}

TEST_F(SimulationManagerTest, TestLiftingLineTwoPanel){
    LiftingSurface ls = LiftingSurface(5,1);
    ls.setAspectRatio( 5.0 );
    ls.setPitch( 10.0 * M_PI / 180.0 );
    ls.getLattice().setHasTrailers(true);   
    ls.getLattice().setTrailerVec( Vec3D(0.0, 1000.0, 0.0) );   
    ls.updateLattice();
    sm->addSurface(&ls);

    LiftingSurface ls2 = LiftingSurface(ls);
    ls2.getLattice().translate( Vec3D(5.0, 0.0, 0.0) );
    sm->addSurface(&ls2);
    
    sm->setGlobalLinearVelocity( Vec3D(0.0, 100.0, 0.0) );
    sm->solve();   
    double gamma[10] = { 35.03083712393401328, 43.13601948056102486, 46.29353105076378938, 47.71129676895412075, 48.28231926132752250, 48.28231926132752960, 47.71129676895412075, 46.29353105076377517, 43.13601948056103197, 35.03083712393401328};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_DOUBLE_EQ( gamma[n], ls.getLattice().getGamma()[i][j] );
        }
    }     
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_DOUBLE_EQ( gamma[n], ls2.getLattice().getGamma()[i][j] );
        }
    }     
}


TEST_F(SimulationManagerTest, TestHighARLiftingLine){
    LiftingSurface ls = LiftingSurface(1,1);
    ls.setAspectRatio( 1.0 );
    ls.setPitch( 5.0 * M_PI / 180.0 );
    ls.getLattice().setHasTrailers(false);   
    ls.getLattice().setTrailerVec( Vec3D(0.0, 10000.0, 0.0) );   
    ls.updateLattice( );
    ls.getLattice().printState();
    sm->addSurface(&ls);
    sm->setGlobalLinearVelocity( Vec3D(0.0, 10.0, 0.0) );
    sm->solve(); 
    printf("Induced Velocity at Gamma: \n");
    ls.getLattice().calcInducedVelocity( ls.getLattice().gammaCenterPoint(0,0) ).printState();
    printf("Cp Normal: \n");
    ls.getLattice().getControlPointNormals()[0][0].printState();
}
