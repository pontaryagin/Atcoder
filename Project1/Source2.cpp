#include "MyHeader.h"
//#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"





// ============================ Header  =================================



int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);


	ll n;
	cin >> n;
	vll h(n), a(n);
	rep(i, 0, n)cin >> h[i];
	rep(i, 0, n)cin >> a[i];

	vll dp(n);
	segment_tree<Monoid::max_t<ll>> seg(n);
	vpll hi(n);
	rep(i, 0, n)hi[i] = pll{ h[i],i };
	sort(all(hi));
	rep(i, 0, n) {
		ll _, j;
		tie(_, j) = hi[i];
		ll m = seg.query(0, j);
		if (m >= 0)
			dp[j] = m + a[j];
		else
			dp[j] = a[j];
		seg.update(j, dp[j]);
	}



	ll res = 0;
	rep(i, 0, n)chmax(res, dp[i]);
	cout << res<<endl;
	return 0;

}
