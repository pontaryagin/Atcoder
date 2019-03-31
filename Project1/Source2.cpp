#include <iostream>
#include <cmath>
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
#include <time.h>
#include <iomanip>
//#include "boost/variant.hpp"

// #include "bits/stdc++.h"

using namespace std;


#define rep(i, N, M) for(ll i=N, i##_len=(M); i<i##_len; ++i)
//#define _overload3(_1, _2, _3, name, ...) name
//#define repi(i,a,b) for(int i=int(a);i<int(b);++i)
//#define _rep(i,n) repi(i,0,n)
//#define rep(...) _overload3(__VA_ARGS__,repi,_rep,)(__VA_ARGS__)

#define rep_skip(i, N, M, ...) for(ll i=N, i##_len=(M); i<i##_len; i+=(skip))
#define rrep(i, N, M)  for(ll i=(M)-1, i##_len=(N-1); i>i##_len; --i)
#define pb push_back

typedef pair<double, double> pd;

typedef long long ll;
typedef pair<ll, ll> pll;

template<int n>
struct tll_impl {
	using type = decltype(tuple_cat(tuple<ll>(), declval<typename tll_impl<n - 1>::type>()));
};
template<>
struct tll_impl<1> {
	using type = tuple<ll>;
};
template<int n> 
using tll = typename tll_impl<n>::type;
// check 
static_assert(is_same<tll<4>, tuple< ll, ll, ll,ll>>::value, "");


template<typename T>
vector<T> make_v(size_t a) { return vector<T>(a); }

template<typename T, typename... Ts>
auto make_v(size_t a, Ts... ts) {
	return vector<decltype(make_v<T>(ts...))>(a, make_v<T>(ts...));
}
// ex:  auto dp =  make_v<ll>(4,5) => vector<vector<ll>> dp(4,vector<ll>(5));

// check if T is vector
template < typename T >
struct is_vector : std::false_type {};

template < typename T >
struct is_vector<vector<T>> : std::true_type {};

template<typename T, typename V>
typename enable_if<!is_vector<T>::value>::type
fill_v(T& t, const V& v) { t = v; }

template<typename T, typename V>
typename enable_if<is_vector<T>::value>::type
fill_v(T& t, const V& v) {
	for (auto &x : t)
		fill_v(x, v);
}
// ex:  fill_v(dp, INF);

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

template<typename T> void chmin(T &a, T b) {
	if (a > b) a = b;
}
template<typename T> void chmax(T &a, T b) {
	if (a < b) a = b;
}

vll seq(ll i, ll j) {
	vll res(j - i);
	rep(k, i, j) res[k] = i + k;
	return res;
}

constexpr ll POW_(ll n, ll m) {
	ll res = 1;
	rep(i, 0, m) {
		res *= n;
	}
	return res;
}

template<ll mod = 0>
constexpr ll POW(ll x, ll n) {
	if (x == 2)
	{
		return (1LL << n) % mod;
	}
	if (n == 0)return 1;
	if (n == 1)return x % mod;
	if (n % 2 == 0)return POW_(POW<mod>(x, n / 2), 2LL) % mod;
	return ((POW_(POW<mod>(x, n / 2), 2LL) % mod)*(x%mod)) % mod;
}
template<>
constexpr ll POW<0>(ll x, ll n) {
	if (x == 2)
	{
		return 1LL << n;
	}
	if (n == 0)return 1;
	if (n == 1)return x;
	if (n % 2 == 0) return POW_(POW(x, n / 2), 2);
	return (POW_(POW(x, n / 2), 2))*x;
}

template<ll mod>
constexpr ll Plus(ll x, ll y) {
	return (x + y) % mod;
}

template<ll mod>
constexpr ll Minus(ll x, ll y) {
	return (x + mod - y) % mod;
}

template<ll mod>
constexpr ll Prod(ll x, ll y) {
	return (x * y) % mod;
}

template<ll mod>
constexpr ll Inv(ll x) {
	assert(x%mod != 0);
	return POW<mod>(x, mod - 2);
}

template<ll mod>
constexpr ll Dev(ll x, ll y) {
	return x * Inv<mod>(y);
}


template<ll bit = 2LL>
ll at_bit(ll n, ll i) {
	return n / POW(bit, i) % bit;
}

template<>
ll at_bit<2>(ll n, ll i) {
	return (n >> i) % 2LL;
}

template<ll bit>
ll get_bit(ll i) {
	return POW(bit, i);
}

template<>
ll get_bit<2>(ll i) {
	return 1LL << i;
}

template<ll bit = 2>
ll get_max_bit(ll n) {
	ll tmp = bit;
	ll at = 0;
	while (tmp <= n) {
		at++;
		tmp *= bit;
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

template<
	typename Inputs,
	typename Functor,
	typename T = typename Inputs::value_type>
	void stable_sort_by(Inputs& inputs, Functor f) {
	std::stable_sort(std::begin(inputs), std::end(inputs),
		[&f](const T& lhs, const T& rhs) { return f(lhs) < f(rhs); });
}

// Functor is expected to be functional<val ,bool>
// This function returns the maximum iterator that stisfies f(*it) == f(* inputs.begin())
template<
	typename Inputs,
	typename Functor,
	typename ValType = typename Inputs::value_type>
	pair<typename Inputs::iterator, typename Inputs::iterator> binary_solve(Inputs& inputs, Functor f)
{

	auto left = inputs.begin();
	auto right = inputs.end();
	auto n = inputs.size();

	auto left_val = f(*left);
	auto right_val = f(*(right - 1));

	// check 
	assert(n >= 2);
	assert(left_val != right_val);

	while (left + 1 != right)
	{
		auto mid = left + (right - left) / 2;
		if (f(*mid) == left_val)
		{
			left = mid;
		}
		else
		{
			right = mid;
		}
	}

	return { left,right };
}

template<typename T>
bool check_binary_solve(T left, T right) {
	return left != right;
}
template<>
bool check_binary_solve<double>(double left, double right) {
	return abs(left - right) > 0.0000000000001;
}

template<typename T, typename Functor>
	pair<T, T > binary_solve(T left, T right, Functor f)
{
	auto left_val = f(left);
	auto right_val = f(right);

	// check 
	assert(left < right);
	if (right_val == left_val)
		throw invalid_argument("right_val == left_val");

	while (check_binary_solve<T>(left, right))
	{
		auto mid = left + (right - left) / 2;
		if (f(mid) == left_val)
		{
			left = mid;
		}
		else
		{
			right = mid;
		}
	}

	return { left,right };
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


constexpr ll MOD = 1000000007;
constexpr ll INF = 1LL << 62;






struct EdgeWithRev { ll to, cap, rev; };
using Weight = ll;
#define RESIDUE(s,t) (capacity[s][t]-flow[s][t])
Weight augment(const Graph &g, const vvll &capacity, vvll &flow,
	const vector<int> &level, vector<bool> &finished, int u, int t, Weight cur) {
	if (u == t || cur == 0) return cur;
	if (finished[u]) return 0;
	finished[u] = true;
	for(auto &e: g[u]) if (level[e.to] > level[u]) {
		Weight f = augment(g, capacity, flow, level, finished,
			e.to, t, min(cur, RESIDUE(u, e.to)));
		if (f > 0) {
			flow[u][e.to] += f; flow[e.to][u] -= f;
			finished[u] = false;
			return f;
		}
	}
	return 0;
}
Weight maximumFlow(const Graph &g, int s, int t) {
	int n = g.size();
	vvll flow(n, vll(n)), capacity(n, vll(n)); // adj. matrix
	rep(u, 0,  n) for(auto e: g[u]) capacity[e.from][e.to] += e.cost;

	Weight total = 0;
	for (bool cont = true; cont; ) {
		cont = false;
		vector<int> level(n, -1); level[s] = 0; // make layered network
		queue<int> Q; Q.push(s);
		for (int d = n; !Q.empty() && level[Q.front()] < d; ) {
			int u = Q.front(); Q.pop();
			if (u == t) d = level[u];
			for(auto &e: g[u]) if (RESIDUE(u, e.to) > 0 && level[e.to] == -1)
				Q.push(e.to), level[e.to] = level[u] + 1;
		}
		vector<bool> finished(n); // make blocking flows
		for (Weight f = 1; f > 0; ) {
			f = augment(g, capacity, flow, level, finished, s, t, INF);
			if (f == 0) break;
			total += f; cont = true;
		}
	}
	return total;
}



int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);

	cout << 234;
	return 0;
}


