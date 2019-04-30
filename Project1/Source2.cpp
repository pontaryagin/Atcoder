
#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
#include "NumberTheory.h"
//#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"



// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	
	ll n;
	cin >> n;
	vll a(n);
	rep(i, 0, n) {
		cin >> a[i];
	}
	sort(all(a));
	reverse(all(a));
	ll one = count(all(a), 1);
	ll two = count(all(a), 2);
	if (one >0  && two>1) {
		a.pop_back();
		a.insert(a.begin() + 1, 1);
	}
	rep(i, 0, n)cout << a[i]<<" ";
	cout << endl;

	return 0;

}
