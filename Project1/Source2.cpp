#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
//#include "NumberTheory.h"
//#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"

// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);

	ll n;
	cin >> n;

	vs M(n);
	rep(i, 0, n) cin >> M[i];

	vector<vector<ll>> N(M.size());
	rep(i, 0, M.size()) {
		rep(j, 0, M[i].size()) {
			N[i].push_back( M[i][j] - '0');

		}
	}

	auto calc = [&](vll & N) {

		vector<ll> a, b;
		if (N.size() == 1) {
			if (N[0] == 5) {
				a.push_back(4);
				b.push_back(1);
			}
			else {
				a.push_back(1);
				b.push_back(N[0] - 1);
			}
		}
		else {
			bool b_on = false;
			if (N[0] == 4) {
				a.push_back(2);
				b.push_back(2);

			}
			else {
				a.push_back(N[0]);
			}
			rep(i, 1, N.size() ) {
				if (N[i] == 5) {
					a.push_back(2);
					b.push_back(3);
				}
				else if(N[i]==0){
					a.push_back(0);
					b.push_back(0);
				}
				else {
					a.push_back(N[i] - 1);
					b.push_back(1);
				}
			}


		}
		return pair<vll, vll>{ a, b };
	};


	rep(i, 0, N.size()) {
		auto res = calc(N[i]);
		cout << "Case #" << i + 1 << ": ";
		rep(j, 0, res.first.size()) {
			cout << res.first[j];
		}
		cout << " ";
		rep(j, 0, res.second.size()) {
			cout << res.second[j];

		}
		cout << endl;
	};

	return 0;
	
}


