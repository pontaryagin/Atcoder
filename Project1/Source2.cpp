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

	ll T;
	cin >> T;
	ll INF_ =numeric_limits<ll>::max() ;
	rep(test, 1, T + 1) {
		[&] {
			ll r, b, c;
			cin >> r >> b >> c;
			vll m(c), s(c), p(c);
			rep(i, 0, c) {
				cin >> m[i] >> s[i] >> p[i];
			}
			ll ng = 0;
			ll ok = INF_;
			while (ng+1 != ok) {
				ll mid = (ng + ok) / 2;
				vll num(c);
				rep(i, 0, c) {
					num[i] = (mid - p[i]) / s[i];
					num[i] = min(num[i], m[i]);
				}
				sort(all(num));
				ll sum = accumulate(num.end() - r, num.end(),0LL);
				if (sum >= b) {
					ok = mid;
				}
				else {
					ng = mid;
				}

			}


			cout << "Case #" << test << ": "<< ok << endl;
		}();
	}

	return 0;

}
