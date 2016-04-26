#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "liftingsurface.h"

class LiftingSurfaceTest: public ::testing::Test {
	public:
	LiftingSurface *ls;
    virtual void SetUp(){
        ls = new LiftingSurface(4,3);
        HorseshoeLattice& hl = ls->getHorseshoeLattice();
        for (int i = 0; i < hl.ni(); i++){
            for (int j = 0; j < hl.nj(); j++){
                hl.getGamma()[i][j] = 1.0;
            }
        }
	}

	virtual void TearDown(){
        delete ls;    
	}
};

TEST_F(LiftingSurfaceTest, TestUpdateLattice){
    ls->setTaperRatio( 0.5 );
    ls->setAspectRatio( 4.0/0.75 ); //Adjusts span
    //vl->snapToAspectTaper( 4.0/0.75, .5);
    ls->updateLattice();    
    HorseshoeLattice& hl = ls->getHorseshoeLattice();
    
    double y[5] =      {0.0,     1.0    ,  2.0    ,  3.0,     4.0};
    double x[4] =      {0.0,    -1.0/3.0, -2.0/3.0, -1.0};
    double xScale[5] = {8.0/8.0, 7.0/8.0, 6.0/8.0, 5.0/8.0, 4.0/8.0};
    double xLe[5] =    {0.0, -(1.0-7.0/8.0)/4.0, -(1.0-6.0/8.0)/4.0, -(1.0-5.0/8.0)/4.0, -(1.0-4.0/8.0)/4.0};

    for (int i = 0; i < hl.ni(); i++){
        for (int j = 0; j < hl.nj(); j++){
           //printf("%i, %i\n",i,j);
           EXPECT_NEAR     (y[i]                      , hl.getEndPoints()[i][j].y, 2E-15 );
           EXPECT_NEAR     (x[j]*xScale[i]+xLe[i]+0.25, hl.getEndPoints()[i][j].x, 2E-15 );
        }
    }
}

TEST_F(LiftingSurfaceTest, SetSetAspectRatio){
    ls->setAspectRatio(8.0);
    EXPECT_DOUBLE_EQ( 8.0, ls->getSpan() );
    
    ls->setAspectRatio( ls->getAspectRatio()/2.0 );
    EXPECT_DOUBLE_EQ( 4.0, ls->getSpan() );
}


