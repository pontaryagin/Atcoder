#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
//#include "NumberTheory.h"
//#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"


// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	vll a(3);
	read_v(a);

	ll H, W, N, sr, sc;
	cin >> H >> W >> N >> sr >> sc;
	sr--; sc--;
	string S, T;
	cin >> S >> T;

	map<char, ll> shift = {
		{'R', 1},{'L',-1},{'U',-1},{'D',1}
	};
	auto drop = [&](ll pos, char dir, char rev) {
		rep(j, 0, N) {
			ll bound = (dir == 'R' || dir == 'L') ? W : H;
			if (S[j] == dir) 
				pos = pos + shift[dir];
			if (pos == -1 || pos == bound) {
				cout << "NO\n"; 
				return true;
			}
			if (T[j] == rev && (pos + shift[rev] >=0 )&& (pos + shift[rev] <= bound -1))
				pos = pos + shift[rev];
		}
		return false;
	};
	if (drop(sc, 'R', 'L') || drop(sc, 'L', 'R') || drop(sr, 'U', 'D') || drop(sr, 'D', 'U')) {

	}
	else cout << "YES\n";
	return 0;

}
