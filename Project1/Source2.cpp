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
	vvll A(n, vll(n));
	vector<tll<3>> data(n*n);
	vector<Edge> edges;
	rep(i, 0, n) {
		rep(j, 0, n) {
			cin >> A[i][j];
		}
	}
	Graph g(A,0);
	auto d = warshall_floyd(g,0);
	ll sum = 0;
	rep(i, 0, n)rep(j, 0, n) {
		if(A[i][j] > d[i][j]) {
			cout << -1 << endl; return 0;
		}
		
	}
	cout << sum<< endl;

	return 0;

}
