#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <math.h>

#include "vec3d.h"

class Vec3DTest: public ::testing::Test {
	public:
    Vec3D *vx, *vy, *vz;    
	virtual void SetUp(){
        vx = new Vec3D(1.0, 0.0, 0.0);
        vy = new Vec3D(0.0, 1.0, 0.0);
        vz = new Vec3D(0.0, 0.0, 1.0);
       
	}

	virtual void TearDown(){
        delete vx; delete vy; delete vz;
	}

};
TEST_F(Vec3DTest, TestCopy){
    Vec3D v1 = Vec3D(1.0, 2.0, 3.0);
    Vec3D v2 = Vec3D(v1);

	EXPECT_DOUBLE_EQ( v2.x, v1.x);
	EXPECT_DOUBLE_EQ( v2.y, v1.y);
	EXPECT_DOUBLE_EQ( v2.z, v1.z);
    
}

TEST_F(Vec3DTest, TestCrossPruduct){
    Vec3D v;
    v  = vx->cross(*vy);
	EXPECT_DOUBLE_EQ( 0.0, v.x);
	EXPECT_DOUBLE_EQ( 0.0, v.y);
	EXPECT_DOUBLE_EQ( 1.0, v.z);
    
    v = vy->cross(*vx);
	EXPECT_DOUBLE_EQ( 0.0, v.x);
	EXPECT_DOUBLE_EQ( 0.0, v.y);
	EXPECT_DOUBLE_EQ( -1.0, v.z);

    v = vy->cross(*vz);
	EXPECT_DOUBLE_EQ( 1.0, v.x);
	EXPECT_DOUBLE_EQ( 0.0, v.y);
	EXPECT_DOUBLE_EQ( 0.0, v.z);
    
    v = vz->cross(*vy);
	EXPECT_DOUBLE_EQ( -1.0, v.x);
	EXPECT_DOUBLE_EQ( 0.0, v.y);
	EXPECT_DOUBLE_EQ( 0.0, v.z);

    v = vz->cross(*vx);
	EXPECT_DOUBLE_EQ( 0.0, v.x);
	EXPECT_DOUBLE_EQ( 1.0, v.y);
	EXPECT_DOUBLE_EQ( 0.0, v.z);

    v = vx->cross(*vz);
	EXPECT_DOUBLE_EQ( 0.0, v.x);
	EXPECT_DOUBLE_EQ( -1.0, v.y);
	EXPECT_DOUBLE_EQ( 0.0, v.z);
    
    Vec3D r1 = Vec3D( 1.0, 1.0, 1.0 );
    Vec3D r2 = Vec3D(-1.0,-1.0, 1.0 );
    Vec3D r3 = r1.cross(r2);
   
    EXPECT_DOUBLE_EQ( 2.0, r3.x );
    EXPECT_DOUBLE_EQ(-2.0, r3.y );
    EXPECT_DOUBLE_EQ( 0.0, r3.z );
    
}

TEST_F(Vec3DTest, TestDotProduct){
    Vec3D v135 = Vec3D(1.0, 3.0, 5.0);
    Vec3D v71113 = Vec3D(7.0, 11.0, 13.0);
    double vDot = v135.dot(v71113);

    EXPECT_DOUBLE_EQ(7.0+33.0+65.0, vDot);
}

TEST_F(Vec3DTest, TestMultiply){
    Vec3D v135 = Vec3D(1.0, 3.0, 5.0);
    double d = 2.0;
    Vec3D dv = d*v135;
    
    EXPECT_DOUBLE_EQ(2.0, dv.x);
    EXPECT_DOUBLE_EQ(6.0, dv.y);
    EXPECT_DOUBLE_EQ(10.0, dv.z);

    Vec3D vd = v135*d;
    
    EXPECT_DOUBLE_EQ(2.0, vd.x);
    EXPECT_DOUBLE_EQ(6.0, vd.y);
    EXPECT_DOUBLE_EQ(10.0, vd.z);
}

TEST_F(Vec3DTest, TestMagnitude){
    Vec3D v135 = Vec3D(1.0, 3.0, 5.0);
    EXPECT_DOUBLE_EQ( 5.916079783099616, v135.magnitude() );
}

TEST_F(Vec3DTest, TestRotationAboutOrigin){
    //Rotation about ones self does nothing
    Vec3D vRot; 
    vRot = vx->rotate(Vec3D(), *vx, M_PI/2.0);
    EXPECT_DOUBLE_EQ( vx->x, vRot.x); 
    EXPECT_DOUBLE_EQ( vx->y, vRot.y); 
    EXPECT_DOUBLE_EQ( vx->z, vRot.z); 
    vRot = vy->rotate(Vec3D(), *vy, M_PI/2.0);
    EXPECT_DOUBLE_EQ( vy->x, vRot.x); 
    EXPECT_DOUBLE_EQ( vy->y, vRot.y); 
    EXPECT_DOUBLE_EQ( vy->z, vRot.z); 
    vRot = vz->rotate(Vec3D(), *vz, M_PI/2.0);
    EXPECT_DOUBLE_EQ( vz->x, vRot.x); 
    EXPECT_DOUBLE_EQ( vz->y, vRot.y); 
    EXPECT_DOUBLE_EQ( vz->z, vRot.z); 
    
    //Rotation of unit vectors about the origin result in unit vectors
    vRot = vx->rotate(Vec3D(), *vy, M_PI/2.0);
    EXPECT_NEAR     ( 0.0, vRot.x, 1E-16);
    EXPECT_NEAR     ( 0.0, vRot.y, 1E-16);
    EXPECT_DOUBLE_EQ(-1.0, vRot.z);
    vRot = vx->rotate(Vec3D(), *vy,-M_PI/2.0);
    EXPECT_NEAR     ( 0.0, vRot.x, 1E-16);
    EXPECT_NEAR     ( 0.0, vRot.y, 1E-16);
    EXPECT_DOUBLE_EQ( 1.0, vRot.z);
    vRot = vy->rotate(Vec3D(), *vz, M_PI/2.0);
    EXPECT_DOUBLE_EQ(-1.0, vRot.x);
    EXPECT_NEAR     ( 0.0, vRot.y, 1E-16);
    EXPECT_NEAR     ( 0.0, vRot.z, 1E-16);
    vRot = vy->rotate(Vec3D(), *vz,-M_PI/2.0);
    EXPECT_DOUBLE_EQ( 1.0, vRot.x);
    EXPECT_NEAR     ( 0.0, vRot.y, 1E-16);
    EXPECT_NEAR     ( 0.0, vRot.z, 1E-16);
    vRot = vz->rotate(Vec3D(), *vx, M_PI/2.0);
    EXPECT_NEAR     ( 0.0, vRot.x, 1E-16);
    EXPECT_DOUBLE_EQ(-1.0, vRot.y);
    EXPECT_NEAR     ( 0.0, vRot.z, 1E-16);
    vRot = vz->rotate(Vec3D(), *vx,-M_PI/2.0);
    EXPECT_NEAR     ( 0.0, vRot.x, 1E-16);
    EXPECT_DOUBLE_EQ( 1.0, vRot.y);
    EXPECT_NEAR     ( 0.0, vRot.z, 1E-16);
}

TEST_F(Vec3DTest, TestRotationAboutPoint){
    //Rotation of unit vectors about non-origin
    Vec3D vRot; 
    vRot = vx->rotate(*vy, *vz, M_PI/2.0);
    EXPECT_NEAR     ( 1.0, vRot.x, 1E-16);
    EXPECT_NEAR     ( 2.0, vRot.y, 1E-16);
    EXPECT_NEAR     ( 0.0, vRot.z, 1E-16);
    vRot = vy->rotate(*vz, *vx, M_PI/2.0);
    EXPECT_NEAR     ( 0.0, vRot.x, 1E-16);
    EXPECT_NEAR     ( 1.0, vRot.y, 1E-16);
    EXPECT_NEAR     ( 2.0, vRot.z, 1E-16);

}
    
TEST_F(Vec3DTest, TestNorm){
    Vec3D v = Vec3D(1.0, 3.0, 5.0);
    Vec3D vNorm = v.norm();
    EXPECT_DOUBLE_EQ(0.1690308509457033, vNorm.x);
    EXPECT_DOUBLE_EQ(0.507092552837110, vNorm.y);
    EXPECT_NEAR     (0.845154254728517, vNorm.z, 1E-15);
}

TEST_F(Vec3DTest, TestSum){
    Vec3D v135 = Vec3D(1.0, 3.0, 5.0);
    double d = 1.0;
    Vec3D dv = d+v135;
    
    EXPECT_DOUBLE_EQ(2.0, dv.x);
    EXPECT_DOUBLE_EQ(4.0, dv.y);
    EXPECT_DOUBLE_EQ(6.0, dv.z);

    Vec3D vd = v135+d;
    
    EXPECT_DOUBLE_EQ(2.0, vd.x);
    EXPECT_DOUBLE_EQ(4.0, vd.y);
    EXPECT_DOUBLE_EQ(6.0, vd.z);

    
    Vec3D v71113 = Vec3D(7.0, 11.0, 13.0);
    Vec3D vv = v135+v71113;
    
    EXPECT_DOUBLE_EQ(8.0, vv.x);
    EXPECT_DOUBLE_EQ(14.0, vv.y);
    EXPECT_DOUBLE_EQ(18.0, vv.z);
}

TEST_F(Vec3DTest, TestDifference){
    Vec3D v135 = Vec3D(1.0, 3.0, 5.0);
    double d = 1.0;
    Vec3D dv = d-v135;
    
    EXPECT_DOUBLE_EQ(0.0, dv.x);
    EXPECT_DOUBLE_EQ(-2.0, dv.y);
    EXPECT_DOUBLE_EQ(-4.0, dv.z);

    Vec3D vd = v135-d;
    
    EXPECT_DOUBLE_EQ(0.0, vd.x);
    EXPECT_DOUBLE_EQ(2.0, vd.y);
    EXPECT_DOUBLE_EQ(4.0, vd.z);

    
    Vec3D v71117 = Vec3D(7.0, 11.0, 17.0);
    Vec3D vv = v135-v71117;
    
    EXPECT_DOUBLE_EQ(-6.0, vv.x);
    EXPECT_DOUBLE_EQ(-8.0, vv.y);
    EXPECT_DOUBLE_EQ(-12.0, vv.z);
}

TEST_F(Vec3DTest, TestAdditionAssignment){
    Vec3D a = Vec3D(1.0, 2.0, 3.0);
    Vec3D b = Vec3D(4.0, 5.0, 6.0);
    a+=b;
    EXPECT_DOUBLE_EQ(5.0, a.x);
    EXPECT_DOUBLE_EQ(7.0, a.y);
    EXPECT_DOUBLE_EQ(9.0, a.z);
    a+=b;
    EXPECT_DOUBLE_EQ(9.0, a.x);
    EXPECT_DOUBLE_EQ(12.0, a.y);
    EXPECT_DOUBLE_EQ(15.0, a.z);
}


