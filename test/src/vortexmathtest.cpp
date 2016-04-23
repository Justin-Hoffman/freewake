#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "vortexmath.h"

class VortexMathTest: public ::testing::Test {
	public:
	virtual void SetUp(){
	}

	virtual void TearDown(){
	}

};

TEST_F(VortexMathTest, TestBiotSavart){
    Vec3D r1 = Vec3D( 1.0, 1.0, 1.0 );
    Vec3D r2 = Vec3D(-1.0,-1.0, 1.0 );
    double rc = 0.0;
    
    Vec3D a = BiotSavart( r1, r2, rc);
	
    EXPECT_DOUBLE_EQ( 0.091888149236965339, a.x ); 
	EXPECT_DOUBLE_EQ(-0.091888149236965339, a.y ); 
	EXPECT_DOUBLE_EQ( 0.000000000000000, a.z );

    r1 = Vec3D(-1.0, 0.0, 0.0 );
    r2 = Vec3D( 0.0, 1.0, 0.0 );
    Vec3D vOut = BiotSavart( r1, r2, 0.0 );
    EXPECT_DOUBLE_EQ( 0.0, vOut.x );
    EXPECT_DOUBLE_EQ( 0.0, vOut.y );
    EXPECT_DOUBLE_EQ( -2.0 / (4.0 * M_PI ) , vOut.z );
     
}
