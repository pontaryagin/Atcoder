#include "pch.h"
#include "../Project1/NumberTheory.h"
#include "../Project1/SegmentTree.h"
#include "../Project1/Graph.h"

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

TEST(Graph, Dijkstra) {
	Graph g(5);
	g.push_undir({ 0,1 });
	g.push_undir({ 0,2, 2 });
	g.push_undir({ 2,3 , 4});
	g.push_undir({ 2,4 });
	vll shortestPathInfo;
	auto depth = g.dijkstra(0, shortestPathInfo);
	auto depth_res = vll{0, 1, 2, 6, 3};
	EXPECT_EQ(depth, depth_res);
	vll path1 = shortest_path_generator(shortestPathInfo, 0, 3);
	auto path1_res = vll{ 0,2,3 };
	EXPECT_EQ(path1, path1_res);
	vll path2 = shortest_path_generator(shortestPathInfo, 2, 3);
	auto path2_res = vll{ 2,3 };
	EXPECT_EQ(path2, path2_res);
}

TEST(Graph, DFSBFS) {
	Graph g(5);
	g.push_undir({ 0,1 });
	g.push_undir({ 1,2 });
	g.push_undir({ 0,3 });
	g.push_undir({ 3,4 });
	vll dfs, bfs;
	g.dfs(0, [&](const Edge & e) {dfs.push_back(e.to); });
	g.bfs(0, [&](const Edge & e) {bfs.push_back(e.to); });
	vll dfs_res = { 1,2,3,4 };
	vll bfs_res = { 1,3,2,4 };
	EXPECT_EQ(dfs, dfs_res);
	EXPECT_EQ(bfs, bfs_res);
}

TEST(Graph, LCA) {

	Graph g(5);
	g.push_undir({ 0,1 });
	g.push_undir({ 0,2 });
	g.push_undir({ 2,3 });
	g.push_undir({ 2,4 });
	LCA lca(g, 0);
	EXPECT_EQ(lca(1, 4), 0);
	EXPECT_EQ(lca(3, 4), 2);
	EXPECT_EQ(lca(2, 4), 2);
	EXPECT_EQ(lca(0, 4), 0);
	LCA lca4(g, 4);
	EXPECT_EQ(lca4(3, 0), 2);
	EXPECT_EQ(lca4(1, 4), 4);

}

TEST(Graph, EulerTour) {

	Graph g(5);
	g.push_undir({ 0,1 });
	g.push_undir({ 0,2 });
	g.push_undir({ 2,3 });
	g.push_undir({ 2,4 });
	// 0
	// \1  \2
	//        \3  \4
	auto tour = euler_tour(g);
	auto res = vll{ 0,1,0,2,3,2,4,2,0 };
	EXPECT_EQ(tour, res);

}
