
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
#include <type_traits>
#include <array>
#include <cassert>
#include <typeinfo>

// #include "bits/stdc++.h"

using namespace std;



#define rep(i, N, M) for(ll i=N, i##_len=(M); i<i##_len; ++i)
#define rep_skip(i, N, M, ...) for(ll i=N, i##_len=(M); i<i##_len; i+=(skip))
#define rrep(i, N, M)  for(ll i=(M)-1, i##_len=(N-1); i>i##_len; --i)
#define pb push_back

typedef pair<double, double> pd;

typedef long long ll;
typedef unsigned long long ull;
typedef pair<ll, ll> pll;
typedef tuple<ll, ll, ll> tll;
typedef tuple<ll, ll, ll, ll> tll4;
typedef tuple<ll, ll, ll, ll, ll> tll5;
typedef tuple<ll, ll, ll, ll, ll, ll> tll6;

typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<pll> vpll;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef vector<string> vs;
template<typename T>
using pq_greater = priority_queue<T, vector<T>, greater<T>>;

#define all(a)  (a).begin(),(a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define vec(a) vector<a>
#define perm(c) sort(all(c));for(bool c##perm=1;c##perm;c##perm=next_permutation(all(c)))

ll POW_(ll n, ll m) {
	ll res = 1;
	rep(i, 0, m) {
		res *= n;
	}
	return res;
}

template<int mod=0>
ll POW(ll x, ll n) {
	if (x == 2)
	{
		return (1LL << n) % mod;
	}
	if (n == 0)return 1;
	if (n == 1)return x % mod;
	if (n % 2 == 0)return POW_(POW<mod>(x, n / 2), 2LL) % mod;
	return ((POW_(POW<mod>(x, n / 2, mod), 2LL) % mod)*(x%mod)) % mod;
}
template<>
ll POW<0>(ll x, ll n) {
	if (x == 2)
	{
		return 1LL << n;
	}
	if (n == 0)return 1;
	if (n == 1)return x ;
	if (n % 2 == 0) return POW_(POW(x, n / 2), 2);
	return (POW_(POW(x, n / 2), 2))*x ;
}

template<ll bit>
ll at_bit(ll n, ll i) {
	return n / POW(bit, i) % bit;
}

template<>
ll at_bit<2>(ll n, ll i) {
	return (n >>i) % 2LL;
}

template<ll bit>
ll get_bit(ll i) {
	return POW(bit,i);
}

template<>
ll get_bit<2>(ll i) {
	return 1LL << i;
}

template<ll bit>
ll get_max_bit(ll n) {
	ll tmp = bit;
	ll at=0;
	while (tmp <= n) {
		at++;
		tmp*= bit;
	}
	return at;
}

template<>
ll get_max_bit<2>(ll n) {
	ll tmp = 2;
	ll at = 0;
	while (tmp <= n) {
		at++;
		tmp <<= 1;
	}
	return at;
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

template<
	typename Inputs,
	typename Functor,
	typename T = typename Inputs::value_type>
void sort_by(Inputs& inputs, Functor f) {
	std::sort(std::begin(inputs), std::end(inputs),
		[&f](const T& lhs, const T& rhs) { return f(lhs) < f(rhs); });
}


template<int I>
vll proj(vpll v) {
	vll res(v.size());
	rep(i, 0, v.size()) {
		if (!I) res[i] = v[i].first;
		else res[i] = v[i].second;
	}
	return res;
}


template<int I, class T>
vll proj(T v) {
	vll res(v.size());
	rep(i, 0, v.size()) {
		res[i] = get<I>(v[i]);
	}
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


ll MOD = 1000000007;
ll INF = 1LL << 60;
ll n;

template <class T> using reversed_priority_queue = priority_queue<T, vector<T>, greater<T> >;

template <typename Monoid>
struct segment_tree 
{

	using underlying_type = typename  Monoid::underlying_type;

	segment_tree(ll a_n) : size_original(a_n)
	{
		vector<underlying_type> initial_value = vector<underlying_type>(a_n, Monoid::unit());
		segment_tree_impl(a_n, initial_value);
	}

	segment_tree(ll a_n, vector<underlying_type>& initial_value) : size_original(a_n)
	{
		segment_tree_impl(a_n, initial_value);
	}

	void update(int i, underlying_type z) { // 0-based
		assert(0 <= i && i < 2 * n - 1);
		a[i + n - 1] = z;
		for (i = (i + n) / 2; i > 0; i /= 2) { // 1-based
			a[i - 1] = Monoid::append(a[2 * i - 1], a[2 * i]);
		}
	}

	underlying_type query(ll l, ll r) { // 0-based, [l, r)
		underlying_type lacc = Monoid::unit(), racc = Monoid::unit();
		assert(l <= r && r <= n);
		l += n; r += n;
		for (; l < r; l /= 2, r /= 2) { // 1-based loop, 2x faster than recursion
			if (l % 2 == 1) lacc = Monoid::append(lacc, a[(l++) - 1]);
			if (r % 2 == 1) racc = Monoid::append(a[(--r) - 1], racc);
		}
		return Monoid::append(lacc, racc);
	}

	ll size() { return size_original; }

private:
	ll size_original;
	ll n;
	vector<underlying_type> a;
	void segment_tree_impl(ll a_n, vector<underlying_type>& initial_value)
	{
		assert(a_n == initial_value.size());
		n = 1; while (n < a_n) n *= 2;
		a.resize(2 * n - 1, Monoid::unit());
		rep(i, 0, initial_value.size()) {
			a[i + (n - 1)] = initial_value[i];
		}
		rrep(i, 0, n - 1) a[i] = Monoid::append(a[2 * i + 1], a[2 * i + 2]); // propagate initial values
	}
};


template <typename T>
struct min_indexed_t {
	typedef pair<T, ll> underlying_type;
	static underlying_type make_indexed(vector<T> v)
	{
		underlying_type w(v.size());
		rep(i, 0, v.size()) {
			w[i] = { v[i],i };
		}
		return w;
	}
	static underlying_type unit() { return make_pair(numeric_limits<T>::max(), -1); }
	static underlying_type append(underlying_type a, underlying_type b) { return min(a, b); }
};
template <typename T>
struct min_t {
	typedef T underlying_type;
	static underlying_type unit() { return numeric_limits<T>::max(); }
	static underlying_type append(underlying_type a, underlying_type b) { return min(a, b); }
};
struct linear_t {
	typedef pd underlying_type;
	static underlying_type unit() { return underlying_type{ 1.,0. }; }
	static underlying_type append(underlying_type a, underlying_type b) {
		return underlying_type{ a.first*b.first, b.first*a.second + b.second };
	}
};




int main() {

	ll N;
	cin >> N;
	vll a(N);
	rep(i, 0, N) {
		cin >> a[i];
	}
	ll neg = count_if(all(a), [](ll x) {return x < 0; });
	ll pos = N - neg;
	ll tmpNeg=0;
	ll tmpPos = neg;
	ll m_ind = 0;
	ll M_ind = 0;

	vll table_m(N+1);
	vll table_M(N+1);
	//fill tabel_m :
	//count of times of <0
	vll b = a;
	rep(i, 0, N + 1) {
		if (i == 0 ) {
			table_m[i] = 0;
		}
		else if (i == 1) {
			table_m[i] = 0;

			if (i > 0 && b[i - 1] >= 0) b[i - 1] *= -2;
		}
		else {
			if (i > 0 && b[i - 1] >= 0) b[i-1]*=-2;
			//calc num i>1
			rrep(j, 1, i) {
				while (b[j] < b[j - 1]) {
					table_m[i]+=2;
					b[j - 1] *=4;
				}
			}
			table_m[i] += table_m[i - 1];
		}
	}

	rrep(i, 0, N + 1) {
		if (i == N) {
			table_M[i] = 0;
		}
		else if (i == N-1) {
			table_M[i] = 0;

			if (i > 0 && b[i - 1] >= 0) b[i - 1] *= -2;
		}
		else {
			if (i > 0 && b[i - 1] >= 0) b[i - 1] *= -2;
			//calc num i>1
			rep(j, i, N-1) {
				while (b[j] > b[j + 1]) {
					table_M[i] += 2;
					b[j +1] *=4;
				}
			}
			table_M[i] += table_M[i + 1];
		}
	}

	vll res(N + 1);
	rep(i, 0, N + 1) {
		ll cnt = 0;
		if (i > 0 && a[i - 1] >= 0) tmpNeg++;
		if (i > 0 && a[i - 1] < 0) tmpPos--;
		cnt += tmpNeg + tmpPos;
		cnt += table_m[i];
		cnt += table_M[i];

		res[i] = cnt;
	}
	cout << *min_element(all(res));

	return 0;
}