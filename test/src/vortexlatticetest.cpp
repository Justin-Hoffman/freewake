#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "vortexlattice.h"

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
