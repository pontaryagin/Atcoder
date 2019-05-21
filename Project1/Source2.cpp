#include "MyHeader.h"
#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
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
	set<pll> me;
	dp[0] = a[0];
	me.insert(pll{ h[0], 0 });
	rep(i, 1, n) {
		auto lb = me.lower_bound(pll{ h[i],i });
		if (lb == me.begin()) {
			dp[i] = a[i];
		}
		else {
			dp[i] = a[i] + dp[(--lb)->second];
		}
	}


	ll res = 0;
	rep(i, 0, n)chmax(res, dp[i]);
	cout << res<<endl;
	return 0;

}
