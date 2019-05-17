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
	cout  << setprecision(12);
	auto x = modint<MOD>(1);
	auto y = modint<>();
	auto z  = y + x;

	ll n, k;
	cin >> n >> k;
	vll a(n);
	read_v(a);
	auto dp = make_v<modint<>>(n, k + 1);
	rep(i, 0, n) {
		if (i == 0)
		{
			rep(j, 0, a[i] + 1)dp[i][j] = 1;
			continue;
		}
		vector<modint<>> S(k + 2);
		rep(j, 0, k + 1)
		{
			S[j + 1] = S[j]; S[j+1] +=dp[i - 1][j];
		}
		rep(j, 0, k+1) {
			if(j ==0) dp[i][0] = 1;
			//dp[i][1] = dp[i - 1][0] + dp[i - 1][1];
			else 
				dp[i][j] = ( S[j + 1] - S[max(0LL, j - a[i])]); //dp[i - 1][max(0LL, j - a[i])] + .. + dp[i - 1][j];
		
		}
	}
	modint<>().test();
	cout << dp[n - 1][k]<<endl ;

	return 0;

}

