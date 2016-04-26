#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "vortexlattice.h"

#include "horseshoelattice.h"
#include "liftingsurface.h"

class VortexLatticeTest: public ::testing::Test {
	public:
	VortexLattice *vl;
    virtual void SetUp(){
        vl = new VortexLattice(6,3);
	}

	virtual void TearDown(){
        delete vl;
	}

};

TEST_F(VortexLatticeTest, TestCopy){
    VortexLattice v2 = VortexLattice(*vl);
    EXPECT_EQ( vl->ni(), v2.ni());    
    EXPECT_EQ( vl->nj(), v2.nj());
}

TEST_F(VortexLatticeTest, TestEndPointsArraySize){
    EXPECT_EQ( 6, vl->endPoints().size() );
    EXPECT_EQ( 3, vl->endPoints()[0].size() );
}

TEST_F(VortexLatticeTest, TestGammaIArraySize){
    EXPECT_EQ( 5, vl->gammaI().size() );
    EXPECT_EQ( 3, vl->gammaI()[0].size() );
}

TEST_F(VortexLatticeTest, TestGammaJArraySize){
    EXPECT_EQ( 6, vl->gammaJ().size() );
    EXPECT_EQ( 2, vl->gammaJ()[0].size() );
}

TEST_F(VortexLatticeTest, TestIJFromN){ 
    //Solutions derived on paper for a 6X3 lattice
    EXPECT_EQ( 0, vl->ijFromN( 0 ).first );
    EXPECT_EQ( 0, vl->ijFromN( 0 ).second );
    
    EXPECT_EQ( 0, vl->ijFromN( 5 ).first );
    EXPECT_EQ( 5, vl->ijFromN( 5 ).second );
    
    EXPECT_EQ( 1, vl->ijFromN( 6 ).first );
    EXPECT_EQ( 0, vl->ijFromN( 6 ).second );
    
    EXPECT_EQ( 1, vl->ijFromN( 7 ).first );
    EXPECT_EQ( 1, vl->ijFromN( 7 ).second );
    
    EXPECT_EQ( 1, vl->ijFromN( 11 ).first );
    EXPECT_EQ( 5, vl->ijFromN( 11 ).second );
    
}

TEST_F(VortexLatticeTest, TestSnapToTrailingEdge){
    LiftingSurface ls = LiftingSurface(4,3);
    ls.setTaperRatio( 0.5 );
    ls.setAspectRatio( 4.0/0.75 ); //Adjusts span
    ls.updateLattice();  
    HorseshoeLattice hl = ls.getHorseshoeLattice();  

    VortexLattice vl = VortexLattice(5,5); //VortexLattice counts edges so we need 5 spanwise instead of 4
    vl.fixToTrailingEdge( hl );
    for (int i = 0; i < vl.ni(); i++){
        EXPECT_DOUBLE_EQ( hl.getEndPoints()[i][3].x, vl.endPoints()[i][0].x );
        EXPECT_DOUBLE_EQ( hl.getEndPoints()[i][3].y, vl.endPoints()[i][0].y );
        EXPECT_DOUBLE_EQ( hl.getEndPoints()[i][3].z, vl.endPoints()[i][0].z );
    } 
}
