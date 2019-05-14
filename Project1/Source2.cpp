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
	cout << fixed << setprecision(12);


	//auto x = prime_factorize(292);
	ll n;
	cin >> n;
	vll a(n);
	rep(i, 0, n)cin >> a[i];
	vll S(n+1);
	rep(i, 0, n) S[i+1] = S[i] ^a[i];
	rep(i, 0, n) cout << S[i+1]<<" ";
	set<ll> all;
	rep(i, 0, n+1)if(S[i]) all.insert(S[i]);
	ll res = 0;
	if (S[n] == 0) {
		vll dp(1LL<<21); // 空間計算量が1<<20程度であることを使っている
		ll sum = 0;
		dp[0]=1;
		rep(i,0,n+1){
			ll x = S[i];
			// search 0x0x0x0;
			if(x==0){
				dp[0]=dp[0]+ sum;
			}else{
				dp[x] = dp[x]+dp[0];
			}
			string s;
		}
		ll u = accumulate(all(dp),0LL);
		cout<< u<< endl;
	}
	else {
		ll x = S[n];
		// 0x0x0x0x;
		ll dp0=1;
		ll dp1=0;

		rep(i, 0, n + 1) {
			if (S[i] != 0) {
				dp1=dp1+dp0;
			}
			else{
				dp0=dp0+dp1;
			}
		}
		cout<< dp1+dp0<<endl;
	}
	return 0;

}

