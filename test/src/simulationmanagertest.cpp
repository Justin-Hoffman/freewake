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
    sm->setReferenceVelocity( 100.0 );
    sm->setReferenceSurface( ReferenceSurface( 10.0, 10.0, 1.0) );
    LiftingSurface ls = LiftingSurface(10,1);
    ls.setAspectRatio( 10.0 );
    ls.setPitch( 0.0 * M_PI / 180.0 );
    ls.getHorseshoeLattice().setHasTrailers(true);   
    ls.getHorseshoeLattice().setTrailerVec( Vec3D(-1000.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -10.0*M_PI/180.0 ) );   
    ls.updateLattice( );
    sm->addSurface(&ls);
    sm->setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -10.0*M_PI/180.0 ) );
    sm->solve();   
    //double gamma[10] = { 35.03083712393401328, 43.13601948056102486, 46.29353105076378938, 47.71129676895412075, 48.28231926132752250, 48.28231926132752960, 47.71129676895412075, 46.29353105076377517, 43.13601948056103197, 35.03083712393401328};
    double gamma[10] = { 35.03083712394531091, 43.13601948057333857, 46.29353105077682073, 47.71129676896759975, 48.28231926134120044, 48.28231926134120044, 47.71129676896760685, 46.29353105077682073, 43.13601948057334567, 35.03083712394531091};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_DOUBLE_EQ( gamma[n], ls.getHorseshoeLattice().getGamma()[i][j] );
        }
    }
}

TEST_F(SimulationManagerTest, TestLiftingLineTwoChord){
    sm->setReferenceVelocity(100.0);
    sm->setReferenceSurface( ReferenceSurface( 10.0, 10.0, 1.0) );
    LiftingSurface ls = LiftingSurface(10,2);
    ls.setAspectRatio( 10.0 );
    ls.setPitch( 10.0 * M_PI / 180.0 );
    ls.getHorseshoeLattice().setHasTrailers(true);   
    ls.getHorseshoeLattice().setTrailerVec( Vec3D(-1000.0, 0.0, 0.0) );   
    ls.updateLattice( );
    sm->addSurface(&ls);
    sm->setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0) );
    sm->solve();   
    //double gamma[10] = { 35.03083712393401328, 43.13601948056102486, 46.29353105076378938, 47.71129676895412075, 48.28231926132752250, 48.28231926132752960, 47.71129676895412075, 46.29353105076377517, 43.13601948056103197, 35.03083712393401328};
    double gamma[10] = { 35.03083712394531091, 43.13601948057333857, 46.29353105077682073, 47.71129676896759975, 48.28231926134120044, 48.28231926134120044, 47.71129676896760685, 46.29353105077682073, 43.13601948057334567, 35.03083712394531091};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        double gammaSum = 0.0;
        n++;
        for(int j = 0; j < ls.nChord(); j++){
            gammaSum += ls.getHorseshoeLattice().getGamma()[i][j];
        }    
        EXPECT_NEAR( gamma[n], gammaSum, gamma[n]*0.01); //Two chordwise station should find net circulation to within 1% of one panel
    }
     
}

TEST_F(SimulationManagerTest, TestBertinCummingsExampleSweptWing){
    sm->setReferenceVelocity(100.0);
    sm->setReferenceSurface( ReferenceSurface( 5.0, 5.0, 1.0) );
    LiftingSurface ls = LiftingSurface(4,1);
    ls.setAspectRatio(2.5);
    ls.setSweep( 45.0 * M_PI / 180.0 );
    ls.getHorseshoeLattice().setHasTrailers(true);   
    ls.getHorseshoeLattice().setTrailerVec( Vec3D(-1000.0, 0.0, 0.0).rotate( Vec3D(0.0,0.0,0.0), Vec3D(0.0,1.0,0.0), -1.0*M_PI/180.0) );   
    ls.updateLattice( );
    sm->addSurface(&ls);

    LiftingSurface ls2 = LiftingSurface( ls );
    ls2.setSweep( -45.0 * M_PI / 180.0 );
    ls2.updateLattice( );
    ls2.getHorseshoeLattice().translate( Vec3D(-2.5, -2.5, 0.0) );
    sm->addSurface(&ls2);
   
    EXPECT_TRUE( ls2.getHorseshoeLattice().hasTrailers() ); 
    sm->setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0).rotate( Vec3D(0,0,0), Vec3D(0.0,1.0,0.0), -1.0*M_PI/180.0) );
    sm->solve();   
    //Taken from 7.49a/b/c/d in sixth edition of Bertin and Cummings - Note the data is given only to 3 sig figs
    double gammaBertinCummings[4] = {2.993628011055486, 3.147147396237819, 3.136181725867653, 2.741417592541655};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_NEAR( gammaBertinCummings[n], ls.getHorseshoeLattice().getGamma()[i][j], 5E-3 );
        }
    }      
    for(int i = 0; i < ls2.nSpan(); i++){
        for(int j = 0; j < ls2.nChord(); j++){
            EXPECT_NEAR( gammaBertinCummings[n], ls2.getHorseshoeLattice().getGamma()[i][j], 5E-3 );
            n--;
        }
    }
    EXPECT_NEAR( 0.0601, sm->netLift()/(1.0/2.0 * 1.0 * 100.0 * 100.0 * 5.0 ), 5E-5);
   
    sm->setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0).rotate( Vec3D(0,0,0), Vec3D(0,1,0), -2.0*M_PI/180.0) );
    sm->solve();   
    EXPECT_NEAR( 0.1202, sm->netLift()/(1.0/2.0 * 1.0 * 100.0 * 100.0 * 5.0 ), 1E-4);
}

TEST_F(SimulationManagerTest, TestLiftingLineTwoPanel){
    sm->setReferenceVelocity(100.0);
    sm->setReferenceSurface( ReferenceSurface( 10.0, 10.0, 1.0) );
    LiftingSurface ls = LiftingSurface(5,1);
    ls.setAspectRatio( 5.0 );
    ls.setPitch( 10.0 * M_PI / 180.0 );
    ls.getHorseshoeLattice().setHasTrailers(true);   
    ls.getHorseshoeLattice().setTrailerVec( Vec3D(-1000.0, 0.0, 0.0) );   
    ls.updateLattice();
    sm->addSurface(&ls);

    LiftingSurface ls2 = LiftingSurface(ls);
    ls2.getHorseshoeLattice().translate( Vec3D(0.0, 5.0, 0.0) );
    sm->addSurface(&ls2);
    
    sm->setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0) );
    sm->solve();   
    //double gamma[10] = { 35.03083712393401328, 43.13601948056102486, 46.29353105076378938, 47.71129676895412075, 48.28231926132752250, 48.28231926132752960, 47.71129676895412075, 46.29353105076377517, 43.13601948056103197, 35.03083712393401328};
    double gamma[10] = { 35.03083712394531091, 43.13601948057333857, 46.29353105077682073, 47.71129676896759975, 48.28231926134120044, 48.28231926134120044, 47.711296768967571, 46.29353105077682073, 43.13601948057334567, 35.03083712394531091};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_DOUBLE_EQ( gamma[n], ls.getHorseshoeLattice().getGamma()[i][j] );
        }
    }     
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_DOUBLE_EQ( gamma[n], ls2.getHorseshoeLattice().getGamma()[i][j] );
        }
    }     
}


TEST_F(SimulationManagerTest, TestHighARLiftingLine){
    sm->setReferenceVelocity( 10.0 );
    sm->setReferenceSurface( ReferenceSurface( 1000.0, 1000.0, 1.0) );
    LiftingSurface ls = LiftingSurface(100,1);
    ls.setAspectRatio( 1000.0 );
    double alpha = 5.0* M_PI / 180.0;
    ls.setPitch( alpha );
    ls.getHorseshoeLattice().setHasTrailers( true );   
    ls.getHorseshoeLattice().setTrailerVec( Vec3D(-10000.0, 0.0, 0.0) );   
    ls.updateLattice( );
    sm->addSurface(&ls);
    sm->setGlobalLinearVelocity( Vec3D(-10.0, 0.0, 0.0) );
    sm->solve(); 
    EXPECT_NEAR( 2.0*M_PI* alpha , sm->netLift()/(1.0/2.0 * 10.0 * 10.0 * 1000.0 ), 5E-3);
}

TEST_F(SimulationManagerTest, TestFreeWake){
    sm->setReferenceVelocity( 100.0 );
    sm->setDt( 0.1 );
    sm->setReferenceSurface( ReferenceSurface( 10.0, 10.0, 1.0) );
    LiftingSurface ls = LiftingSurface(10,1,10,2);
    ls.setFreeWake( true );
    ls.setAspectRatio( 10.0 );
    ls.setPitch( 0.0 * M_PI / 180.0 );
    ls.getHorseshoeLattice().setHasTrailers(false);   
    ls.updateLattice( );
    sm->addSurface(&ls);
    sm->setGlobalLinearVelocity( Vec3D(-100.0, 0.0, 0.0).rotate(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), -10.0*M_PI/180.0 ) );
    double Cl = -10.0;
    int nStep = 0;
    while ( nStep < 100 ){
        sm->step();
        nStep++; 
    }
    //double gamma[10] = { 35.03083712393401328, 43.13601948056102486, 46.29353105076378938, 47.71129676895412075, 48.28231926132752250, 48.28231926132752960, 47.71129676895412075, 46.29353105076377517, 43.13601948056103197, 35.03083712393401328};
    double gamma[10] = { 35.03083712394531091, 43.13601948057333857, 46.29353105077682073, 47.71129676896759975, 48.28231926134120044, 48.28231926134120044, 47.71129676896760685, 46.29353105077682073, 43.13601948057334567, 35.03083712394531091};
    int n = -1;
    for(int i = 0; i < ls.nSpan(); i++){
        for(int j = 0; j < ls.nChord(); j++){
            n++;
            EXPECT_NEAR( gamma[n], ls.getHorseshoeLattice().getGamma()[i][j], 2E-2);
        }
    }
}
