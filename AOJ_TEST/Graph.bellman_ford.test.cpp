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

#define PROBLEM "http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_1_B"

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll n, m, st; cin >> n >> m >> st;
	Digraph g(n);
	rep(i, 0, m) {
		Edge e;
		cin >> e.from >> e.to >> e.cost;
		g.push(e);
	}
	auto dist = g.bellman_ford(st);
	rep(i, 0, n) {
		if (dist[i] == -INF) {
			cout << "NEGATIVE CYCLE" << endl;
			return 0;
		}
	}
	rep(i, 0, n) {
		if (dist[i] == INF) {
			cout << "INF" << endl;
		}
		else {
			cout << dist[i] << endl;
		}
	}
	return 0;

}
