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
	vll S(n+1),zero(n+1);
	rep(i, 0, n) { S[i + 1] = S[i] ^ a[i]; zero[i + 1] = zero[i] + (S[i] == 0?1:0); }
	//rep(i, 0, n) cout << S[i+1]<<" ";
	//set<ll> all;
	//rep(i, 0, n+1)if(S[i]) all.insert(S[i]);
	ll res = 0;
	if (S[n] == 0) {
		vll dp0(1LL<<21), dp1(1LL << 21), cnt(1LL<<21), ind(1LL<<21); // 空間計算量が1<<20程度であることを使っている
		ll sum = 0;
		rep(i,1,n+1){
			ll x = S[i];
			// search 0x0x0x0;
			if(x==0){
				//rep(i, 0, dp0.size()) {
					//dp0[i] = dp0[i] + dp1[i]; dp0[i] %= MOD;
				//}
				sum++;
			}else{
				dp0[x] = dp0[x] + ( dp1[x] * (sum-cnt[x]) );
				dp0[x] %= MOD;
				ind[x] = i;
				cnt[x] = sum;
				dp1[x] = dp1[x]+ dp0[x]+1;
				dp1[x] %= MOD;
			}
		}

		ll u=0;
		rep(i, 0, dp1.size()) { u += dp1[i];u %= MOD; }
		ll cn = (count(all(S), 0) - 2);
		u += POW(2, cn, MOD); //POW<MOD>(2LL, cn) ;
		u %= MOD;
		//write_v(S);
		//rep(i, 0, S.size())cout << dp[S[i]];
		cout<< u << endl;
	}
	else {
		ll x = S[n];
		// 0x0x0x0x;
		ll dp0=1;
		ll dp1=0;
		//if(n>100)abort();
		rep(i, 0, n) {
			if (S[i] != 0) {
				if (S[i] == x)
				{
					dp1 = dp1 + dp0; dp1 %= MOD;
				}
					
			}
			else{
				{dp0 = dp0 + dp1; dp0 %= MOD; }
			}
		}
		cout<< dp0<<endl;
	}
	return 0;

}

