#include "MyHeader.h"
#include "Graph.h"
#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"


// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);
	
	ll n, k; cin >> n >> k;
	Graph g(n);
	vll a(n - 1), b(n - 1);
	rep(i, 0, n - 1) {
		cin >> a[i] >> b[i];
		g.push({ a[i],b[i] });
		g.push({ b[i],a[i] });
	}
	GraphDFS dfs(g);
	ll res = 0;
	dfs([&](Edge& e) {
		if (g.out(e.to).size > k - 1) {
			auto sz = g.out;

		}
		}
		, 0);



	return 0;

}