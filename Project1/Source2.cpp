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

	vector<bitset<21>> as(n);
	rep(i, 0, n)rep(j, 0, n) {
		as[i][j] = a[i][j];
	}
	vector<vector<int>> dp(n, vector<int>(1LL << 22));
	rep(j, 0, 1 << n) {
		ll i = popcnt(j)-1;
		auto poss = as[i] & bitset<21>(j);
		if (i == 0) {
			dp[0][j] = poss.count();
			continue;
		}
		if (poss == 0)continue;
		rep(k, 0, n) {
			if (poss[k]) dp[i][j] += dp[i - 1][j ^ (1LL << k)];
			if (dp[i][j] >= MOD) dp[i][j] -= MOD;
		}
	}
	cout << dp[n - 1][(1LL << n) - 1] << endl;
	return 0;

}
