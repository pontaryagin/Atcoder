#include "pch.h"
#include "../Project1/NumberTheory.h"
#include "../Project1/SegmentTree.h"


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

TEST(SegmentTree, LazySegmentTree) {

	LazySegmentTree<M::min_indexed_t<>> seg(10);
	seg.update(0, pll{ 20,0 });
	seg.update(1, pll{ 11,1 });
	seg.update(3, 4, pll{ 3,0 });
	EXPECT_EQ(seg.query(0, 1).first, 20);
	EXPECT_EQ(seg.query(1, 4).first, 3);
	EXPECT_EQ(seg.query(0, 2).first, 11);
	EXPECT_EQ(seg.query(0, 0).first, numeric_limits<ll>::max());

}

TEST(SegmentTree, segment_tree) {

	segment_tree<M::min_indexed_t<>> seg(10);
	seg.update(0, pll{ 20,0 });
	seg.update(1, pll{ 11,1 });
	seg.update(3, pll{ 3,0 });
	EXPECT_EQ(seg.query(0, 1).first, 20);
	EXPECT_EQ(seg.query(1, 4).first, 3);
	EXPECT_EQ(seg.query(0, 2).first, 11);
	EXPECT_EQ(seg.query(0, 0).first, numeric_limits<ll>::max());

}