#include <gtest/gtest.h>
#include <gmock/gmock.h>

class VanillaTest: public ::testing::Test {
	public:
	virtual void SetUp(){
	}

	virtual void TearDown(){
	}

};

TEST_F(VanillaTest, TestTrue){
	EXPECT_TRUE(true); // Continue testing if this fails
	ASSERT_TRUE(true); // Hard stop if this fails
}
