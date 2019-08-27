#include "pch.h"
#include "../Project1/NumberTheory.h"
#include "../Project1/SegmentTree.h"
#include "../Project1/Graph.h"
#include "../Project1/Text.h"

TEST(MyHeader, inv_map) {
	vll x = { 4, 5, 10 };
	auto mp = inv_map(x);
	EXPECT_EQ(mp[4], 0);
	EXPECT_EQ(mp[10], 2);
	vector<double> y = { 0.1, 10., 100., 1000. };
	auto mp2 = inv_map(y);
	EXPECT_EQ(mp2[100.], 2);
}

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

TEST(NumberTheory, runtime_modint) {
	modint<-1>::set_modulo(7);
	auto x = modint<-1>(5);
	auto y = modint<-1>(3);
	modint<-2>::set_modulo(5);
	auto u = modint<-2>(5);
	auto v = modint<-2>(3);
	EXPECT_EQ(x + y, 1);
	EXPECT_EQ(u + v, 3);
	auto w = modint<-3>(11, 5);
	EXPECT_EQ(w, 1);
	auto a = modint<-4>(6);
	a.reset_modulo(5);
	EXPECT_EQ(a, 1);

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
	EXPECT_EQ(seg.query(3).first, 3);

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

TEST(Graph, EdgeItr) {
	Graph g(3);
	g.push({ 0,1, 1});
	g.push({ 0,2, 2 });
	g.push({ 1,2 , 3 });
	ll i = 0;
	for (auto& e : g) {
		EXPECT_EQ(e.cost, ++i);
	}
	EXPECT_EQ(i, 3);
}

TEST(Graph, Bellman_Ford) {
	// normal case
	Graph g(5);
	g.push_undir({ 0,1 });
	g.push_undir({ 0,2, 2 });
	g.push_undir({ 2,3 , 4 });
	g.push_undir({ 2,4 });
	vll shortestPathInfo;
	auto depth = g.bellman_ford(0, shortestPathInfo);
	auto depth_res = vll{ 0, 1, 2, 6, 3 };
	EXPECT_EQ(depth, depth_res);
	vll path1 = shortest_path_generator(shortestPathInfo, 0, 3);
	auto path1_res = vll{ 0,2,3 };
	EXPECT_EQ(path1, path1_res);
	vll path2 = shortest_path_generator(shortestPathInfo, 2, 3);
	auto path2_res = vll{ 2,3 };
	EXPECT_EQ(path2, path2_res);
	// negative loop
	Graph g2(5);
	g2.push({ 0,1, -1 });
	g2.push({ 1,1, -100 });
	g2.push({ 0,3,-1 });
	g2.push({ 2,3,-1 });
	g2.push({ 2,2,-100 });
	g2.push({ 1,4,-1 });
	auto dist = g2.bellman_ford(0,-INF);
	vll res = { 0,-INF,INF,-1,-INF };
	EXPECT_EQ(dist, res);
}

TEST(Graph, DFSBFS) {
	Graph g(5);
	g.push_undir({ 0,1 });
	g.push_undir({ 1,2 });
	g.push_undir({ 0,3 });
	g.push_undir({ 3,4 });
	vll dfs, bfs, dfs_node;
	g.dfs(0, [&](const Edge & e) {dfs.push_back(e.to); });
	g.bfs(0, [&](const Edge & e) {bfs.push_back(e.to); });
	g.dfs_node(0, [&](ll node) {dfs_node.push_back(node); });
	vll dfs_res = { 1,2,3,4 };
	vll bfs_res = { 1,3,2,4 };
	vll dfs_node_res = { 0,1,2,3,4 };
	EXPECT_EQ(dfs, dfs_res);
	EXPECT_EQ(bfs, bfs_res);
	EXPECT_EQ(dfs_node, dfs_node_res);
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
	//       \3  \4
	auto tour = g.euler_tour(0);
	auto res = vll{ 0,1,0,2,3,2,4,2,0 };
	EXPECT_EQ(tour, res);

}

TEST(Graph, Kruskal) {

	Graph g(5);
	g.push_undir({ 0,1,1 });
	g.push_undir({ 0,2,3 });
	g.push_undir({ 0,3,4 });
	g.push_undir({ 1,3,2 });
	g.push_undir({ 2,3,1 });
	// 0 - 1
	// | \ |
	// \2- \3
	Graph mst = g.kruskal();
	set<Edge> edgeMST, edgeRes;
	rep(i, 0, mst.edges.size())
		edgeMST.insert(mst[i]);
	edgeRes.insert({ 0,1,1 });
	edgeRes.insert({ 1,0,1 });
	edgeRes.insert({ 1,3,2 });
	edgeRes.insert({ 3,1,2 });
	edgeRes.insert({ 2,3,1 });
	edgeRes.insert({ 3,2,1});
	EXPECT_EQ(edgeMST, edgeRes);
}

TEST(Graph, is_bipartite) {
	Graph g(5);
	g.push_undir({ 0,1 });
	g.push_undir({ 0,2 });
	g.push_undir({ 0,3 });
	g.push_undir({ 1,3 });
	EXPECT_EQ(g.is_bipartite(), false);
	Graph g2(5);
	g2.push_undir({ 0,1 });
	g2.push_undir({ 1,2 });
	g2.push_undir({ 2,3 });
	g2.push_undir({ 3,0 });
	EXPECT_EQ(g2.is_bipartite(), true);
}

TEST(Graph, max_flow){
	Graph g(6);
	g.push({0,1});
	g.push({0,2});
	g.push({1,2});
	g.push({1,3});
	g.push({1,4});
	g.push({2,3});
	g.push({3,5});
	g.push({4,5});
	FordFulkerson ff(g);
	EXPECT_EQ(ff.max_flow(0, 5), 2);
}

TEST(Text, RollingHash) {
	auto rh = RollingHash<>("abcdefg");
	auto rh2 = RollingHash<>("xbc");
	EXPECT_EQ(rh(1, 3), rh2(1, 3));
}

TEST(Text, kmp_search) {
	EXPECT_EQ(kmp_search("aaaabbaabcb", "aabc"),6);
	EXPECT_EQ(kmp_search("aaaaaa", "aaa"), 0);
	EXPECT_EQ(kmp_search("aaa", "aaaa"),3);
}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
