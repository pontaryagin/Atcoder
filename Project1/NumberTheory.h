#pragma once
#include "MyHeader.h"




ll div_ferm(ll a, ll  b, ll mod) {
	return (a* POW(b, mod - 2, mod)) % mod;
}


// === Modint ===

template <std::uint_fast64_t Modulus = 1000000007> 
class modint 
{
	using u64 = std::uint_fast64_t;

public:
	u64 a;

	constexpr modint(const u64 x = 0) noexcept : a(x % Modulus) {}
	constexpr modint(const modint& rhs) noexcept {
		this->a = rhs.value();
	}
	constexpr modint &operator=(const modint &rhs) noexcept {
		this->a = rhs.value();
		return *this;
	}
	constexpr u64 value() const noexcept { return a; }
	constexpr modint operator+(const modint rhs) const noexcept {
		return modint(*this) += rhs;
	}
	constexpr modint operator-(const modint rhs) const noexcept {
		return modint(*this) -= rhs;
	}
	constexpr modint operator*(const modint rhs) const noexcept {
		return modint(*this) *= rhs;
	}
	constexpr modint operator/(const modint rhs) const noexcept {
		return modint(*this) /= rhs;
	}
	modint &operator+=(const modint rhs) noexcept {
		a += rhs.a;
		if (a >= Modulus) {
			a -= Modulus;
		}
		return *this;
	}
	modint &operator-=(const modint rhs) noexcept {
		if (a < rhs.a) {
			a += Modulus;
		}
		a -= rhs.a;
		return *this;
	}
	modint &operator*=(const modint rhs) noexcept {
		a = a * rhs.a % Modulus;
		return *this;
	}
	modint &operator/=(modint rhs) noexcept {
		u64 exp = Modulus - 2;
		while (exp) {
			if (exp % 2) {
				*this *= rhs;
			}
			rhs *= rhs;
			exp /= 2;
		}
		return *this;
	}
	modint &operator++() noexcept {
		return *this += modint(1);
	}
	modint &operator++(int) noexcept {
		auto t = *this;
		*this += modint(1);
		return t;
	}
	modint &operator--() noexcept {
		return *this -= modint(1);
	}
	modint &operator--(int) noexcept {
		auto t = *this;
		*this -= modint(1);
		return t;
	}
	constexpr bool operator==(const modint rhs) const noexcept { return a == rhs.value(); }
	constexpr bool operator!=(const modint rhs)const  noexcept { return a != rhs.value(); }
	constexpr bool operator <(const modint rhs)const  noexcept { return a < rhs.value(); }
	// should be moved to Google Test
	//constexpr void test(){
	//	constexpr auto x = modint<5>(3);
	//	constexpr auto y = modint<5>(4);
	//	constexpr auto z = modint<5>(2);
	//	static_assert(x + y == z, "");
	//	static_assert(x != y, "");
	//	static_assert(++x == y && x++ != y && x == --y && x != y--, "");
	//	static_assert(x + 6 == y, "");
	//	static_assert(x / 2 == y, "");
	//	
	//}
};

template<uint_fast64_t Modulus>
ostream& operator <<(ostream &o, const modint<Modulus> &t) {
	o << t.value();
	return o;
}
template<uint_fast64_t Modulus>
istream& operator >>(istream &in, modint<Modulus> &t) {
	uint_fast64_t x;
	in >> x;
	t = modint<Modulus>(x);
	return in;
}
template<uint_fast64_t Modulus>
modint<Modulus> POW(modint<Modulus> x, ll n) {
	return modint<Modulus>(POW(x.value(), n, Modulus));
}

// === Mll ===

template<typename T, T MOD = 1000000007>
struct Mint {
	T v;
	Mint() :v(0) {}
	Mint(signed v) :v(v) {}
	Mint(long long t) { v = t % MOD; if (v < 0) v += MOD; }

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

	static Mint pow(Mint v, long long k) {
		Mint res(1), tmp(v);
		while (k) {
			if (k & 1) res *= tmp;
			tmp *= tmp;
			k >>= 1;
		}
		return res;
	}

	// find x s.t. a^x = b
	static T log(Mint a, Mint b) {
		const T sq = 40000;
		unordered_map<T, T> dp;
		dp.reserve(sq);
		Mint res(1);
		for (ll r = 0; r < sq; r++) {
			if (!dp.count(res)) dp[res] = r;
			res *= a;
		}
		Mint p = pow(a.inv(), sq);
		res = b;
		for (ll q = 0; q <= MOD / sq + 1; q++) {
			if (dp.count(res)) {
				T idx = q * sq + dp[res];
				if (idx > 0) return idx;
			}
			res *= p;
		}
		return T(-1);
	}

	static vector<Mint> fact, finv, invs;

	static void init(ll n) {
		if (n + 1 <= (signed)fact.size()) return;
		fact.assign(n + 1, 1);
		finv.assign(n + 1, 1);
		invs.assign(n + 1, 1);

		for (ll i = 1; i <= n; i++) fact[i] = fact[i - 1] * Mint(i);
		finv[n] = Mint(1) / fact[n];
		for (ll i = n; i >= 1; i--) finv[i - 1] = finv[i] * Mint(i);
		for (ll i = 1; i <= n; i++) invs[i] = finv[i] * fact[i - 1];
	}

	static Mint comb(long long n, ll k) {
		Mint res(1);
		for (ll i = 0; i < k; i++) {
			res *= Mint(n - i);
			res /= Mint(i + 1);
		}
		return res;
	}

	static Mint C(ll n, ll k) {
		if (n < k || k < 0) return Mint(0);
		init(n);
		return fact[n] * finv[n - k] * finv[k];
	}

	static Mint P(ll n, ll k) {
		if (n < k || k < 0) return Mint(0);
		init(n);
		return fact[n] * finv[n - k];
	}

	static Mint H(ll n, ll k) {
		if (n < 0 || k < 0) return Mint(0);
		if (!n && !k) return Mint(1);
		init(n + k - 1);
		return C(n + k - 1, k);
	}

	static Mint S(ll n, ll k) {
		Mint res;
		init(k);
		for (ll i = 1; i <= k; i++) {
			Mint tmp = C(k, i)*Mint(i).pow(n);
			if ((k - i) & 1) res -= tmp;
			else res += tmp;
		}
		return res *= finv[k];
	}

	static vector<vector<Mint> > D(ll n, ll m) {
		vector<vector<Mint> > dp(n + 1, vector<Mint>(m + 1, 0));
		dp[0][0] = Mint(1);
		for (ll i = 0; i <= n; i++) {
			for (ll j = 1; j <= m; j++) {
				if (i - j >= 0) dp[i][j] = dp[i][j - 1] + dp[i - j][j];
				else dp[i][j] = dp[i][j - 1];
			}
		}
		return dp;
	}

	static Mint B(ll n, ll k) {
		Mint res;
		for (ll j = 1; j <= k; j++) res += S(n, j);
		return res;
	}

	static Mint montmort(ll n) {
		Mint res;
		init(n);
		for (ll k = 2; k <= n; k++) {
			if (k & 1) res -= finv[k];
			else res += finv[k];
		}
		return res *= fact[n];
	}

	static Mint LagrangePolynomial(vector<Mint> &y, Mint t) {
		ll n = y.size() - 1;
		if (t.v <= n) return y[t.v];
		init(n + 1);
		Mint num(1);
		for (ll i = 0; i <= n; i++) num *= t - Mint(i);
		Mint res;
		for (ll i = 0; i <= n; i++) {
			Mint tmp = y[i] * num / (t - Mint(i))*finv[i] * finv[n - i];
			if ((n - i) & 1) res -= tmp;
			else res += tmp;
		}
		return res;
	}
};


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
		:MOD(MOD), MAX(max(MAX, 2LL)), fac(vll(MAX + 1)), finv(vll(MAX + 1)), inv(vll(MAX + 1)) {
		fac[0] = fac[1] = 1;
		finv[0] = finv[1] = 1;
		inv[1] = 1;
		pre_process(2LL, MAX + 1);
	}

	ll get(ll n, ll k) {
		if (MAX < n)
			pre_process(MAX + 1, n + 1);

		if (n < k)return 0;
		if (n < 0 || k < 0)return 0;
		return fac[n] * (finv[k] * finv[n - k] % MOD) % MOD;
	}
private:
	void pre_process(ll m, ll n) {
		if (MAX < n) {
			fac.resize(n); inv.resize(n); finv.resize(n);
		}
		rep(i, m, n) {
			fac[i] = fac[i - 1] * i % MOD;
			inv[i] = MOD - inv[MOD%i] * (MOD / i) % MOD;
			finv[i] = finv[i - 1] * inv[i] % MOD;
		}
	}
};


ll choose(int n, int r) { // O(r) for small n
	ll acc = 1;
	rep(i, 0, r) acc = acc * (n - i) / (i + 1);
	return acc;
}

ll gcd(ll a, ll b) {
	if (a%b == 0) return b;
	else return gcd(b, a%b);
}




vll getDivisors(ll n) {
	vll res;
	ll i = 1;

	for (; i*i < n; i++) {
		if (n%i == 0) {
			res.push_back(i);
			res.push_back(n / i);
		}
	}
	if (i*i == n)res.push_back(i);
	sort(res.begin(), res.end());
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




