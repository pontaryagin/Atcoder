#pragma once

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
#define rep_skip(i, N, M, ...) for(ll i=N, i##_len=(M); i<i##_len; i+=(skip))
#define rrep(i, N, M)  for(ll i=(M)-1, i##_len=(N-1); i>i##_len; --i)
#define pb push_back

constexpr ll MOD = 1000000007;
constexpr ll INF = 1LL << 62;

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


template<int n, typename T>
struct vec_t_impl {
	using type = vector<typename vec_t_impl<n - 1, T>::type>;
};
template<typename T>
struct vec_t_impl<1, T> {
	using type = vector<T>;
};
template<int n, typename T>
using vec_t = typename vec_t_impl<n, T>::type;
// check 
static_assert(is_same<vec_t<3, ll>, vector<vector<vector<ll>>>>::value, "");

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

// ============================ Header  =================================




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


