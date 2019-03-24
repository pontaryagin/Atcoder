#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
#include <cstdio>
#include <string>
#include <bitset>
#include <list>
#include <set>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <functional>
#include <queue>
#include <regex>
#include <cassert>
#include <map>
// #include "bits/stdc++.h"

using namespace std;
//debug

#define rep(i, N, M) for (ll i = N; i < M; ++i)
#define rrep(i, N, M) for (ll i = N; i > M; --i)
#define pb push_back


typedef long long ll;
typedef unsigned long long ull;
typedef pair<ll, ll> pll;
typedef pair<int, int> pii;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef vector<string> vs;
typedef priority_queue<ll> pqll;
typedef priority_queue<pll, vector<pll>> pqpll;
typedef priority_queue<ll, vll, greater<ll>> pqll_greater;
typedef priority_queue<pll, vector<pll>, greater<pll>> pqpll_greater;


#define all(a)  (a).begin(),(a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define vec(a) vector<a>
#define _CRT_SECURE_NO_WARNINGS



ll POW(ll n, ll m) {
	ll res = 1;
	rep(i, 0, m) {
		res *= n;
	}
	return res;
}

ll bin_power(ll x, ll y, ll mod) {
	if (y == 0)return 1;
	if (y == 1)return x % mod;
	if (y % 2 == 0)return POW(bin_power(x, y / 2, mod), 2LL) % mod;
	return ((POW(bin_power(x, y / 2, mod), 2LL) % mod)*(x%mod)) % mod;
}
ll div_ferm(ll a, ll  b, ll mod) {
	return (a* bin_power(b, mod - 2, mod)) % mod;
}


vector<pll >prime_factorize(ll n) {
	vector<pll> res;
	for (ll p = 2; p*p <= n; ++p) {
		if (n%p != 0)continue;
		ll num = 0;
		while (n%p == 0) { ++num; n /= p; }
		res.push_back({ p,num });
	}
	if (n != 1) res.push_back(make_pair(n, 1));
	return res;
}

class Combination {
	// this calculates combination (nCk).
	// Constructor runs in O(MAX).
	// get(n,k) returns nCk in O(1).

	ll MAX, MOD;
	vll fac;
	vll finv;
	vll inv;
public:
	Combination(ll MAX = 210000, ll MOD = 1000000007)
		:MOD(MOD), MAX(MAX), fac(vll(MAX)), finv(vll(MAX)), inv(vll(MAX)) {
		fac[0] = fac[1] = 1;
		finv[0] = finv[1] = 1;
		inv[1] = 1;
		rep(i, 2, MAX) {
			fac[i] = fac[i - 1] * i % MOD;
			inv[i] = MOD - inv[MOD%i] * (MOD / i) % MOD;
			finv[i] = finv[i - 1] * inv[i] % MOD;
		}
	}

	ll get(ll n, ll k) {
		if (n < k)return 0;
		if (n < 0 || k < 0)return 0;
		return fac[n] * (finv[k] * finv[n - k] % MOD) % MOD;
	}
};


ll choose(int n, int r) { // O(r) for small n
	ll acc = 1;
	rep(i, 0, r) acc = acc * (n - i) / (i + 1);
	return acc;
}

vll getDivisors(ll n) {
	// O(sqrt(n))
	vll res;
	ll i = 1;
	for (; i*i < n; i++) {
		if (n%i == 0) {
			res.push_back(i);
			res.push_back(n / i);
		}
	}
	if (i*i == n)res.push_back(i);
	return res;
}

vll getDivisors(ll n, ll m) {
	// O(sqrt(min(n,m)))
	if (n > m) swap(n, m);
	vll res;
	ll i = 1;

	for (; i*i < n; i++) {
		if (n%i == 0) {
			if (m%i == 0) res.push_back(i);
			if (m % (n / i) == 0) res.push_back(n / i);
		}
	}
	if (i*i == n) if (m%i == 0) res.push_back(i);
	sort(res.begin(), res.end());
	return res;
}

ll gcd(ll a, ll b) {
	if (a%b == 0) return b;
	else return gcd(b, a%b);
}

template<typename T, T MOD = 1000000007>
struct Mint {
	T v;
	Mint() :v(0) {}
	Mint(signed v) :v(v) {}
	Mint(long long t) { v = t % MOD; if (v < 0) v += MOD; }

	Mint pow(long long k) {
		Mint res(1), tmp(v);
		while (k) {
			if (k & 1) res *= tmp;
			tmp *= tmp;
			k >>= 1;
		}
		return res;
	}

	Mint inv() { return pow(MOD - 2); }

	Mint& operator+=(Mint a) { v += a.v; if (v >= MOD)v -= MOD; return *this; }
	Mint& operator-=(Mint a) { v += MOD - a.v; if (v >= MOD)v -= MOD; return *this; }
	Mint& operator*=(Mint a) { v = 1LL * v*a.v%MOD; return *this; }
	Mint& operator/=(Mint a) { return (*this) *= a.inv(); }

	Mint operator+(Mint a) const { return Mint(v) += a; };
	Mint operator-(Mint a) const { return Mint(v) -= a; };
	Mint operator*(Mint a) const { return Mint(v) *= a; };
	Mint operator/(Mint a) const { return Mint(v) /= a; };

	Mint operator-() { return v ? MOD - v : v; }

	bool operator==(const Mint a)const { return v == a.v; }
	bool operator!=(const Mint a)const { return v != a.v; }
	bool operator <(const Mint a)const { return v < a.v; }

	// find x s.t. a^x = b
	static T log(Mint a, Mint b) {
		const T sq = 40000;
		unordered_map<T, T> dp;
		dp.reserve(sq);
		Mint res(1);
		for (Int r = 0; r < sq; r++) {
			if (!dp.count(res)) dp[res] = r;
			res *= a;
		}
		Mint p = pow(a.inv(), sq);
		res = b;
		for (Int q = 0; q <= MOD / sq + 1; q++) {
			if (dp.count(res)) {
				T idx = q * sq + dp[res];
				if (idx > 0) return idx;
			}
			res *= p;
		}
		return T(-1);
	}

	static vector<Mint> fact, finv, invs;

	static void init(Int n) {
		if (n + 1 <= (signed)fact.size()) return;
		fact.assign(n + 1, 1);
		finv.assign(n + 1, 1);
		invs.assign(n + 1, 1);

		for (Int i = 1; i <= n; i++) fact[i] = fact[i - 1] * Mint(i);
		finv[n] = Mint(1) / fact[n];
		for (Int i = n; i >= 1; i--) finv[i - 1] = finv[i] * Mint(i);
		for (Int i = 1; i <= n; i++) invs[i] = finv[i] * fact[i - 1];
	}

	static Mint comb(long long n, Int k) {
		Mint res(1);
		for (Int i = 0; i < k; i++) {
			res *= Mint(n - i);
			res /= Mint(i + 1);
		}
		return res;
	}

	static Mint C(Int n, Int k) {
		if (n < k || k < 0) return Mint(0);
		init(n);
		return fact[n] * finv[n - k] * finv[k];
	}

	static Mint P(Int n, Int k) {
		if (n < k || k < 0) return Mint(0);
		init(n);
		return fact[n] * finv[n - k];
	}

	static Mint H(Int n, Int k) {
		if (n < 0 || k < 0) return Mint(0);
		if (!n && !k) return Mint(1);
		init(n + k - 1);
		return C(n + k - 1, k);
	}

	static Mint S(Int n, Int k) {
		Mint res;
		init(k);
		for (Int i = 1; i <= k; i++) {
			Mint tmp = C(k, i)*Mint(i).pow(n);
			if ((k - i) & 1) res -= tmp;
			else res += tmp;
		}
		return res *= finv[k];
	}

	static vector<vector<Mint> > D(Int n, Int m) {
		vector<vector<Mint> > dp(n + 1, vector<Mint>(m + 1, 0));
		dp[0][0] = Mint(1);
		for (Int i = 0; i <= n; i++) {
			for (Int j = 1; j <= m; j++) {
				if (i - j >= 0) dp[i][j] = dp[i][j - 1] + dp[i - j][j];
				else dp[i][j] = dp[i][j - 1];
			}
		}
		return dp;
	}

	static Mint B(Int n, Int k) {
		Mint res;
		for (Int j = 1; j <= k; j++) res += S(n, j);
		return res;
	}

	static Mint montmort(Int n) {
		Mint res;
		init(n);
		for (Int k = 2; k <= n; k++) {
			if (k & 1) res -= finv[k];
			else res += finv[k];
		}
		return res *= fact[n];
	}

	static Mint LagrangePolynomial(vector<Mint> &y, Mint t) {
		Int n = y.size() - 1;
		if (t.v <= n) return y[t.v];
		init(n + 1);
		Mint num(1);
		for (Int i = 0; i <= n; i++) num *= t - Mint(i);
		Mint res;
		for (Int i = 0; i <= n; i++) {
			Mint tmp = y[i] * num / (t - Mint(i))*finv[i] * finv[n - i];
			if ((n - i) & 1) res -= tmp;
			else res += tmp;
		}
		return res;
	}
};

