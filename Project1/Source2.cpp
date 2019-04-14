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

	mt19937 rng(5);


	
	ll T;
	cin >> T;
	vll vR(T), vC(T);


	rep(k, 0, T) {
		cin >> vR[k] >> vC[k];

	}
	rep(k, 0, T) {
		ll R = vR[k], C = vC[k];
		auto check = [&](ll x, ll y) {
			ll r = x / C;
			ll rr = y / C;
			ll c = x % C;
			ll cc = y % C;

			return !(r == rr || c == cc || r - c == rr - cc || r + c == rr + cc);
		};
		vll res;
		set<ll> visited;
		vvll data(R*C);
		rep(i, 0, R * C)rep(j, 0, R * C) {
			if(check(i,j))
				data[i].push_back(j);

		}
		rep(i, 0, R* C) {
			random_shuffle(data[i].begin(), data[i].end());
		}
		bool dfs_finished = false;
		function<bool(ll)> dfs = [&](ll st) {
			if (res.size() == R * C) {
				return true;
			}

			for(ll i: data[st]) {

				if (visited.find(i) == visited.end()) {
					visited.insert(i);
					res.push_back(i);
					if (dfs(i))
						return true;
				}
			}
			visited.erase(st);
			res.pop_back();
			return false;
		};

		cout << "Case #" << k + 1 << ": ";

		[&]() {
			rep(i, 0, R * C) {
				visited.insert(i);
				res.push_back(i);
				dfs(i);
				if (dfs(i)) {
					cout << "POSSIBLE" << endl;
					rep(i, 0, res.size()) {
						cout << res[i] / C + 1 << " " << res[i] % C + 1 << endl;
					}
					return;
				}
			}
			cout << "IMPOSSIBLE" << endl;
		}();
		
	}
	

	return 0;
	
}
