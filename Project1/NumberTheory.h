#pragma once
#include "MyHeader.h"




ll div_ferm(ll a, ll  b, ll mod) {
	return (a* POW(b, mod - 2, mod)) % mod;
}


// === Modint ===

template <std::uint_fast64_t Modulus = MOD> 
class modint 
{
	using u64 = std::uint_fast64_t;

public:
	u64 a;

	constexpr modint(const u64 x = 0) noexcept : a(x % Modulus) {}
	modint inv() { return pow(Modulus - 2); }
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
	constexpr modint operator-() { return a ? Modulus - a : a; }
	constexpr bool operator==(const modint rhs) const noexcept { return a == rhs.value(); }
	constexpr bool operator!=(const modint rhs)const  noexcept { return a != rhs.value(); }
	constexpr bool operator <(const modint rhs)const  noexcept { return a < rhs.value(); }
	static constexpr modint zero() { return modint(0); }
	static constexpr modint unit() { return modint(1); }

	modint pow(long long k) const {
		modint v = *this;
		modint res(1), tmp(v);
		while (k) {
			if (k & 1) res *= tmp;
			tmp *= tmp;
			k >>= 1;
		}
		return res;
	}
	
	u64 log(modint b) {
		modint a = *this;
		const u64 sq = 40000;
		map<modint, u64> dp;
		//dp.reserve(sq);
		modint res(1);
		for (ll r = 0; r < sq; r++) {
			if (!dp.count(res)) dp[res] = r;
			res *= a;
		}
		modint p = a.inv().pow(sq);
		res = b;
		for (ll q = 0; q <= Modulus / sq + 1; q++) {
			if (dp.count(res)) {
				u64 idx = q * sq + dp[res];
				if (idx > 0) return idx;
			}
			res *= p;
		}
		return INF;
	}

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



class Combination {
	// this calculates combination (nCk).
	// Constructor runs in O(MAX).
	// get(n,k) returns nCk in O(1).

	ll N_MAX, mod;
	vll fac;
	vll finv;
	vll inv;
public:
	Combination(ll mod = MOD, ll N_MAX = 210000)
		:mod(mod), N_MAX(max(N_MAX, 2LL)), fac(vll(N_MAX + 1)), finv(vll(N_MAX + 1)), inv(vll(N_MAX + 1)) {
		fac[0] = fac[1] = 1;
		finv[0] = finv[1] = 1;
		inv[1] = 1;
		pre_process(2LL, N_MAX + 1);
	}

	ll operator()(ll n, ll k) {
		if (N_MAX < n)
			pre_process(N_MAX + 1, n + 1);

		if (n < k)return 0;
		if (n < 0 || k < 0)return 0;
		return fac[n] * (finv[k] * finv[n - k] % mod) % mod;
	}
private:
	void pre_process(ll m, ll n) {
		if (N_MAX < n) {
			fac.resize(n); inv.resize(n); finv.resize(n);
		}
		rep(i, m, n) {
			fac[i] = fac[i - 1] * i % mod;
			inv[i] = mod - inv[mod%i] * (mod / i) % mod;
			finv[i] = finv[i - 1] * inv[i] % mod;
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




