#include "pch.h"
#include "../Project1/NumberTheory.h"

TEST(NumberTheory,modint) {

	auto x = modint<5>(3);
	auto y = modint<5>(4);
	auto z = modint<5>(2);
	EXPECT_EQ(x + y, z);
	EXPECT_NE(x, y);
	EXPECT_EQ(x + 6 , y);
	EXPECT_EQ(x / 2, y);
	auto w = x;
	EXPECT_TRUE(++w == y);
	w = x;
	EXPECT_TRUE(w++ != y);
	w = y;
	EXPECT_TRUE( x == --w );
	w = y;
	EXPECT_TRUE(x != w--);
	EXPECT_EQ(-x, 2);
	EXPECT_EQ(x.pow( 2), 4);
	EXPECT_EQ(x.pow( 4*100000), 1);
	EXPECT_EQ(x.log( 4), 2);



}

TEST(NumberTheory, Combination){

	Combination cmb;
	EXPECT_EQ(cmb(4, 2), 6);


}