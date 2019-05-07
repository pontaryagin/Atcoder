#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
//#include "NumberTheory.h"
#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"

// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);


	ll n;
	cin >> n;
	vpll ab(n - 1);
	rep(i, 0, n - 1) {
		cin >> ab[i].first >> ab[i].second;
		ab[i].first--; ab[i].second--;
		ab.push_back({ ab[i].second, ab[i].first });
	}

	Graph g(n, ab, Graph::nondir);
	ll dia = g.diameter();
	if (dia % 3 == 1) {
		cout << "Second" << endl;
	}
	else {
		cout << "First" << endl;
	}
	
	return 0;

}
