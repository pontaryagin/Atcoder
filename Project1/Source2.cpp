#include "MyHeader.h"
#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"





// ============================ Header  =================================



int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);


	ll n;
	cin >> n;
	Graph g(n);

	rep(i, 0, n - 1) {
		Edge e;
		cin >> e.from >> e.to >> e.cost;
		e.from--; e.to--;
		g.push_undir(e);
	}

	vll color(n, -1);
	color[0] = 0;
	auto func = [&](Edge& e) {
		ll x = e.to;
		ll dist = e.cost;
		ll col = color[e.from];
		if (dist % 2 == 0) {
			color[x] = col;
		}
		else {
			col ^= 1;
			color[x] = col;

		}
	};
	
	GraphDFS dfs(g);
	dfs(func);
	rep(i, 0, n) {
		cout << color[i] << endl;
	}

	return 0;

}
