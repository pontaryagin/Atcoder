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
		for (auto x : all) {
			// search 0x0x0x0;

		}
	}
	else {
		ll x = S[n];
		// 0x0x0x0x;
		rep(i, 0, n + 1) {
			if (S[i] != 0) {

			}
		}
	}

	return 0;

}

