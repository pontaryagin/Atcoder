#include "MyHeader.h"
//#include "Graph.h"
#include "NumberTheory.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"





// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << setprecision(12);
	auto x = modint<MOD>(1);
	auto y = modint<>();
	auto z = y + x;

	ll n;
	cin >> n;
	vvll a(n, vll(n));
	read_v(a);

	vector<ll> as(n);
	rep(i, 0, n)rep(j, 0, n) {
		as[i] |= 1<< a[i][j];
	}
	vector<vector<Mint<>>> dp(n, vector<Mint<>>(1LL << 22));
	rep(j, 0, 1 << n) {
		ll i = popcnt(j)-1;
		auto poss = as[i] & (j);
		if (i == 0) {
			dp[0][j] = popcnt(poss);
			continue;
		}
		if (poss == 0)continue;
		rep(k, 0, n) {
			if (poss & (1<<k)) dp[i][j] += dp[i - 1][j ^ (1LL << k)];
		}
	}
	cout << dp[n - 1][(1LL << n) - 1].v << endl;
	return 0;

}
