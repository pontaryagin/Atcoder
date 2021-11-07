#include "gtest/gtest.h"
#include "NumberTheory.h"
#include "SegmentTree.h"
#include "Graph.h"
#include "Text.h"
#include "Polynomial.h"
#include "Matrix.h"
#include "RecurrenceRelation.h"
#include "D2.h"

TEST(MyHeader, inv_map) {
    vll x = { 4, 5, 10 };
    auto mp = inv_map(x);
    EXPECT_EQ(mp[4], 0);
    EXPECT_EQ(mp[10], 2);
    vector<double> y = { 0.1, 10., 100., 1000. };
    auto mp2 = inv_map(y);
    EXPECT_EQ(mp2[100.], 2);
    map<ll, ll> mp_o;
    mp_o[2] = 100; mp_o[3] = 121;
    auto mp_i = inv_map(mp_o);
    EXPECT_EQ(mp_i[100], 2);
    EXPECT_EQ(mp_i[121], 3);

}

TEST(MyHeader, complex) {
    complex<ll> c1{ 1,2 };
    complex<ll> c2{ 2,3 };
    EXPECT_EQ(c1 + c2, complex<ll>(3, 5));
    EXPECT_EQ(c1.imag(), 2);
    c2.imag(4);
    EXPECT_EQ(c2, complex<ll>(2, 4));
}
TEST(MyHeader, exist) {
    vector<ll> v{ 0,1,2 };
    set<ll> s{ 0,1,2 };
    map<ll, ll> m{ {0,10},{1,10},{2,20} };
    EXPECT_EQ(exist(v, 2ll), true);
    EXPECT_EQ(exist(s, 2ll), true);
    EXPECT_EQ(exist(m, 2ll), true);
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
    using mod = modint<7>;
    Combination<mod> cmb;
    EXPECT_EQ(cmb(4, 2), 6);
    EXPECT_EQ(cmb.Fac(3), 6);
    EXPECT_EQ(cmb.FacInv(3), mod(1) / 6);
    EXPECT_EQ(cmb.P(3, 3), 6);
    EXPECT_EQ(cmb.H(3, 2), cmb(3 + 2 - 1, 2));

    Combination<double> cmb2;
    EXPECT_EQ(cmb2.Fac(4), 4. * 3. * 2. * 1.);
}

TEST(NumberTheory, divisors) {
    auto divs = divisor(12);
    vll res{ 1,2,3,4,6,12 };
    EXPECT_EQ(divs, res);
    EXPECT_EQ(divisor(4), (vll{ 1,2,4 }));
    auto divs2 = divisor(12, 15);
    vll res2{ 1,3 };
    EXPECT_EQ(divs2, res2);
}

TEST(NumberTheory, prime_factorize) {
    auto primes = prime_factorize(12);
    vpll res{ {2,2},{3,1} };
    EXPECT_EQ(primes, res);
}

TEST(NumberTheory, eulears_phi) {
    EXPECT_EQ(eulers_phi(6), 2);
    EXPECT_EQ(eulers_phi(1000000), 400000);
}

TEST(SegmentTree, LazySegmentTree) {

    LazySegmentTree<M::min_indexed_t<>> seg(10);
    seg.update(0, pll{ 20,0 });
    seg.update(1, pll{ 11,1 });
    seg.update(3, 4, pll{ 3,0 });
    EXPECT_EQ(seg.query(0, 1).first, 20);
    EXPECT_EQ(seg.query(0, 1).second, 0);
    EXPECT_EQ(seg.query(1, 4).first, 3);
    EXPECT_EQ(seg.query(1, 4).second, 0);
    EXPECT_EQ(seg.query(0, 2).first, 11);
    EXPECT_EQ(seg.query(0, 2).second, 1);
    EXPECT_EQ(seg.query(0, 0).first, numeric_limits<ll>::max());

}
TEST(SegmentTree, segment_tree) {
    segment_tree < M::min_t<>> seg(10);
    seg.update(0, 20);
    seg.update(1, 11);
    seg.update(3, 3);
    EXPECT_EQ(seg.query(0, 1), 20);
    EXPECT_EQ(seg.query(1, 4), 3);
    EXPECT_EQ(seg.query(0, 2), 11);
    EXPECT_EQ(seg.query(0, 0), numeric_limits<ll>::max());
    EXPECT_EQ(seg.query(3), 3);

}


TEST(Graph, Dijkstra) {
    Graph g(5);
    g.push({ 0,1 });
    g.push({ 0,2, 2 });
    g.push({ 2,3 , 4});
    g.push({ 2,4 });
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
    vll path3 = shortest_path_generator(shortestPathInfo, 0, 0);
    auto path3_res = vll{ 0 };
    EXPECT_EQ(path3, path3_res);
}

TEST(Graph, EdgeItr) {
    Digraph g(3);
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
    g.push({ 0,1 });
    g.push({ 0,2, 2 });
    g.push({ 2,3 , 4 });
    g.push({ 2,4 });
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
    Digraph g2(5);
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
    g.push({ 0,1 });
    g.push({ 1,2 });
    g.push({ 0,3 });
    g.push({ 3,4 });
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
    g.push({ 0,1 });
    g.push({ 0,2 });
    g.push({ 2,3 });
    g.push({ 2,4 });
    LCA lca(g, 0);
    EXPECT_EQ(lca(1, 4), 0);
    EXPECT_EQ(lca(3, 4), 2);
    EXPECT_EQ(lca(2, 4), 2);
    EXPECT_EQ(lca(0, 4), 0);
    LCA lca4(g, 4);
    EXPECT_EQ(lca4(3, 0), 2);
    EXPECT_EQ(lca4(1, 4), 4);

}

TEST(Graph, EulerTourEdge) {

    Graph g(5);
    g.push({ 0,1 });
    g.push({ 0,2 });
    g.push({ 2,3 });
    g.push({ 2,4 });
    // 0
    // \1  \2
    //       \3  \4
    auto tour = g.euler_tour_edge(0);
    auto res = vll{ 0,1,0,2,3,2,4,2,0 };
    EXPECT_EQ(tour, res);

}

TEST(Graph, EulerTourNode) {

    Graph g(5);
    g.push({ 0,1 });
    g.push({ 0,2 });
    g.push({ 2,3 });
    g.push({ 2,4 });
    auto tour = g.euler_tour_node(0);
    auto res = vvll{ {0, 9}, {1, 2}, {3, 8}, {4, 5}, {6, 7} };
    EXPECT_EQ(tour, res);

}

TEST(Graph, Kruskal) {

    Graph g(5);
    g.push({ 0,1,1 });
    g.push({ 0,2,3 });
    g.push({ 0,3,4 });
    g.push({ 1,3,2 });
    g.push({ 2,3,1 });
    // 0 - 1
    // | \ |
    // \2- \3
    auto mst = g.kruskal<GraphDir::undir>();
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
    g.push({ 0,1 });
    g.push({ 0,2 });
    g.push({ 0,3 });
    g.push({ 1,3 });
    EXPECT_EQ(g.is_bipartite(), false);
    Graph g2(5);
    g2.push({ 0,1 });
    g2.push({ 1,2 });
    g2.push({ 2,3 });
    g2.push({ 3,0 });
    EXPECT_EQ(g2.is_bipartite(), true);
}

TEST(Graph, max_flow){
    Digraph g(6);
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

TEST(Graph, diameter) {
    Graph_Base<GraphDir::undir ,double> g(5);
    g.push({ 0,1,1. });
    g.push({ 1,2,2. });
    g.push({ 2,4,2.5 });
    g.push({ 1,3, .5 });
    EXPECT_EQ(g.diameter(), 5.5);
}

TEST(Graph, max_length) {
    Digraph g(5);
    g.push({ 1,0 });
    g.push({ 1,2 });
    g.push({ 2,3 });
    g.push({ 4,1 });
    EXPECT_EQ(g.max_length(), 3);
}

TEST(Graph, acyclic) {
    Digraph g(5);
    g.push({ 1,0 });
    g.push({ 1,2 });
    g.push({ 2,3 });
    g.push({ 4,1 });
    EXPECT_EQ(g.acyclic(), true);
    g.push({ 4,3 });
    EXPECT_EQ(g.acyclic(), true);
    g.push({ 3,1 });
    EXPECT_EQ(g.acyclic(), false);
    vll loop;
    g.acyclic(loop);
    vll res{ 1,2,3 };
    EXPECT_EQ(loop, res);
}

TEST(Graph, acyclic_undir) {
    Graph g(5);
    g.push({ 1,0 });
    g.push({ 1,2 });
    g.push({ 2,3 });
    g.push({ 4,1 });
    EXPECT_EQ(g.acyclic(), true);
    g.push({ 4,3 });
    EXPECT_EQ(g.acyclic(), false);
    vll res{ 1,2,3,4 };
    vll loop;
    g.acyclic(loop);
    EXPECT_EQ(loop, res);
    g.push({ 3,1 });
    EXPECT_EQ(g.acyclic(), false);
}

TEST(Graph, warshall_floyd) {
    Graph g(4);
    g.push({ 0,1 });
    g.push({ 0,2 });
    auto d = g.warshall_floyd();
    EXPECT_EQ(d[0][1], 1);
    EXPECT_EQ(d[1][2], 2);
    EXPECT_EQ(d[0][3], INF);
    Digraph g2(4);
    g2.push({ 0,1 });
    g2.push({ 0,2 });
    d = g2.warshall_floyd();
    EXPECT_EQ(d[0][1], 1);
    EXPECT_EQ(d[1][2], INF);
    EXPECT_EQ(d[0][3], INF);
}


TEST(Graph, count_loops) {
    Graph g(20);
    g.push({ 0,1 });
    g.push({ 0,2 });
    g.push({ 1,2 });

    g.push({ 4,5 });

    g.push({ 6, 7 });
    g.push({ 6, 8 });
    g.push({ 7, 8 });
    g.push({ 8, 9 });
    g.push({ 9, 10 });
    g.push({ 9, 11 });
    g.push({ 10, 11 });
    
    EXPECT_EQ(g.count_loops(0), 1);
    EXPECT_EQ(g.count_loops(1), 1);

    EXPECT_EQ(g.count_loops(3), 0);

    EXPECT_EQ(g.count_loops(4), 0);

    EXPECT_EQ(g.count_loops(6), 2);
    EXPECT_EQ(g.count_loops(8), 2);
    EXPECT_EQ(g.count_loops(9), 2);
    EXPECT_EQ(g.count_loops(11), 2);
}


TEST(UnionFind, weight) {
    UnionFind uf(4);
    uf.unite(0, 3, 2);
    uf.unite(1, 3, 4);
    EXPECT_EQ(uf.is_same(0, 1), true);
    EXPECT_EQ(uf.is_same(3, 1), true);
    EXPECT_EQ(uf.is_same(1, 2), false);
    EXPECT_EQ(uf.weight_diff(3, 1), -4);
    EXPECT_EQ(uf.weight_diff(0, 3), 2);
    EXPECT_EQ(uf.weight_diff(0, 1), -2);
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

TEST(Text, run_length) {
    auto compressed = run_length(string("aaabcdde"));
    vector<pair<char,ll>> res = { {'a',3}, {'b',1}, {'c',1},{'d',2},{'e',1} };
    EXPECT_EQ(compressed, res);
}

TEST(Polynomial, SparsePolinomial) {
    using Pol = SparsePolynomial<>;
    Pol p({ {1,2}, {0, 1} }); // 1 + 2 * x
    Pol q({ {2,3}, {1,1}, {0, 2} }); // 2 + x + 3*x^2
    auto prod = p * q; // 2 + 5 * x + 5 * x^2 + 6 * x^3
    auto prod_res = Pol({ {0,2}, {1,5}, {2,5}, {3,6} });
    EXPECT_EQ(prod.coeff, prod_res.coeff);
    EXPECT_EQ((q * p).coeff, prod_res.coeff);
    vll prod_each = { prod[0],  prod[1], prod[2], prod[3] };
    vll prod_each_res = vll{ 2, 5, 5, 6 };
    EXPECT_EQ(prod_each, prod_each_res);

    auto sum = p + q; // 3 + 3 * x + 3 * x^2
    auto sum_res = Pol({ {0,3}, {1,3}, {2,3} });
    EXPECT_EQ(sum.coeff, sum_res.coeff);
    EXPECT_EQ((q + p).coeff, sum_res.coeff);

}

TEST(Polynomial, SparsePolunomialWithDim) {
    using Mod = modint<2>;
    using Pol = SparsePolynomial<Mod, 1>;
    Pol p({ {0, 1}, {1,2} }); // 1 + 2 * x
    Pol q({ {0, 2}, {1,1}  }); // 2 + x + 3*x^2
    auto prod = p * q; // 2 + 5 * x + 5 * x^2 + 6 * x^3
    auto prod_res = Pol({ {1,1} });
    EXPECT_EQ(prod.coeff, prod_res.coeff);
    EXPECT_EQ((q*p).coeff, prod_res.coeff);
    auto sum = p + q; // 1 + 3 * x + 3 * x^2
    auto sum_res = Pol({ {0,1}, {1,1}});
    EXPECT_EQ(sum.coeff, sum_res.coeff);
    EXPECT_EQ((q+p).coeff, sum_res.coeff);
    EXPECT_EQ(p << 1, Pol({ {1,1} }));
}

TEST(Matrix, POW) {
    Matrix<ll> m{ 3 }, n{ 3 };
    rep(i, 0, 3) {
        rep(j, 0, 3) {
            m[i][j] = i + j;
            n[i][j] = i * j;
        }
    }
    EXPECT_EQ(m.POW(3), m* m* m);
    EXPECT_EQ(m.POW(4), m* m* m* m);
    EXPECT_EQ(n.POW(3), n* n* n);
    EXPECT_EQ(n.POW(8), n* n* n* n * n* n* n* n);
    EXPECT_EQ(m.POW(1), m);
    EXPECT_EQ(m.POW(0), Matrix<ll>::I(3));

}

TEST(RecurrenceRelation, solve_recurrence_relation) {
    vll p = { 3,2 };
    vll a = { 1, 0 };
    vll a_ans = { 1, 0 };
    a_ans.resize(100);
    rep(i, 2, 100) {
        a_ans[i] = p[0] * a_ans[i - 2] + p[1] * a_ans[i - 1];
    }
    vll b(100);
    rep(i, 0, 100) {
        b[i] = solve_recurrence_relation(a, p, i);
    }
    EXPECT_EQ(a_ans, b);
}

TEST(Polynomial, base_operation) {
    Polynomial<ll> p({1,2});
    Polynomial<ll> q({3,2,1});
    EXPECT_EQ(p * q, Polynomial<ll>({3,8,5,2}));
    EXPECT_EQ(p + q, Polynomial<ll>({ 4,4,1 }));
    EXPECT_EQ((p >> 1).a, Polynomial<ll>({ 2 }).a);
    EXPECT_EQ((p << 1).a, Polynomial<ll>({ 0,1,2 }).a);

}

TEST(Polynomial, div) {
    Polynomial<ll> p({ 1 });
    Polynomial<ll> q({ 1,-1 });
    auto d = p.div(q);
    rep(i, 0, 100) {
        EXPECT_EQ(d(i), 1);
    }
}


TEST(Polynomial, div2) {
    Polynomial<double> p({ 1,2 });
    Polynomial<double> q({ 3,2,1 });
    auto d1 = p.div(q);
    auto d2 = q.div(p);
    vector<double> ans1 = { 1. / 3., 4. / 9., -11. / 27., 10./ 81., 13./243., -56. / 729. };
    vector<double> ans2 = { 3., -4., 9., -18., 36., -72. };
    rep(i, 0, ans1.size()) {
        EXPECT_NEAR(d1(i), ans1[i], 1e-10);
    }
    rep(i, 0, ans2.size()) {
        EXPECT_EQ(d2(i), ans2[i]);
    }
}

TEST(D2, main) {

    ll h = 2, w = 3;
    D2 d2(h, w);
    EXPECT_EQ(d2.next(0, D2::L), nullopt);
    EXPECT_EQ(d2.next(0, D2::D), nullopt);
    EXPECT_EQ(d2.next(0, D2::U), 3);
    EXPECT_EQ(d2.next(0, D2::R), 1);
    EXPECT_EQ(d2(5), pll(1,2));
    EXPECT_EQ(d2(1,2), 5);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
