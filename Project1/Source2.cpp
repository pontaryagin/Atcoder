
#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
//#include "NumberTheory.h"
//#include "Graph.h"
#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"



// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	




	ll T;
	cin >> T;
	rep(test, 1, T + 1) {
		[&] {
			ll n, k;
			cin >> n >> k;
			vll c(n), d(n);
			segment_tree<Monoid::max_t<ll> >  seg(n);

			rep(i, 0, n) {
				cin >> c[i]; 
				seg.update(i, c[i]);
			}
			rep(i, 0, n) { 
				cin >> d[i]; 
				seg.update(i, d[i]);
			}
			ll cnt = 0;

			rep(i, 0, n) {
				rep(j, i, n) {
					if (abs(seg.query(i,j) - seg.query(i,j)) <= k) {
						cnt++;
						
					}
				}
			}

			cout << "Case #" << test << ": " << cnt << endl;
		}();
	}



	return 0;

}
