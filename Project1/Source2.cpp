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
	ll T;
	cin >> T;
	vll x;
	cout << x[T];
}

//int main() {
//	cin.tie(0);
//	ios::sync_with_stdio(false);
//
//	ll T;
//	cin >> T;
//	rep(test, 1, T + 1) {
//		[&] {
//			ll r, b, c;
//			cin >> r >> b >> c;
//			vll m(c), s(c), p(c);
//			rep(i, 0, c) {
//				cin >> m[i] >> s[i] >> p[i];
//			}
//			ll ng = 0;
//			ll ok = numeric_limits<ll>::max();
//			while (ng + 1 != ok) {
//				ll mid = (ng + ok) / 2;
//				vll num;
//				rep(i, 0, c) {
//					ll x  = (mid - p[i]) / s[i];
//					num.push_back(min(x, m[i]));
//				}
//				sort(all(num));
//				reverse(all(num));
//				ll sum = 0;
//				rep(i, 0, min(ll(num.size()), r))sum += num[i];
//				if (sum >= b) {
//					ok = mid;
//				}
//				else {
//					ng = mid;
//				}
//
//			}
//
//
//			cout << "Case #" << test << ": " << ok << endl;
//		}();
//	}
//
//	return 0;
//
//}
