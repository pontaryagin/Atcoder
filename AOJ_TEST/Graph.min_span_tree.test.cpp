#include "MyHeader.h"
#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "Text.h"
//#include "2D.h"

// ============================ Header  =================================

#define PROBLEM "http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_2_A"

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll n, m; cin >> n >> m;
	Graph g(n);
	rep(i, 0, m) {
		ll s, t, w; cin >> s >> t >> w;
		g.push({ s,t,w });
	}
	auto tree = g.kruskal<GraphDir::dir>();
	ll res = 0;
	rep(i, 0, tree.edges.size()) {
		res += tree[i].cost;
	}
	cout << res << endl;
	return 0;

}