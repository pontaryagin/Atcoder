
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



namespace {

#define rep(i, N, M) for(int i=N, i##_len=(M); i<i##_len; ++i)
#define rrep(i, N, M)  for(int i=(M)-1, i##_len=(N-1); i>i##_len; --i)
#define pb push_back


	typedef long long ll;
	typedef unsigned long long ull;
	typedef pair<ll, ll> pll;
	typedef tuple<ll, ll, ll> tll;
	typedef tuple<ll, ll, ll, ll> tll4;
	typedef vector<ll> vll;
	typedef vector<vll> vvll;
	typedef vector<pll> vpll;
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
#define FIR first
#define SEC second
#define perm(c) sort(all(c));for(bool c##perm=1;c##perm;c##perm=next_permutation(all(c)))


	template<class T, class S>
	T atbit(T n, S i) {
		return (n >> i) % i;
	}

	template<class T>
	T getbit(T i) {
		return 1LL << i;
	}
	template<class T>
	T POW(T n, T m) {
		T res = 1;
		rep(i, 0, m) {
			res *= n;
		}
		return res;
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
		vll res;
		ll i = 1;

		for (; i*i < n; i++) {
			if (n%i == 0) {
				if (m%i == 0) res.push_back(i);
				if (m%i == 0) res.push_back(n / i);
			}
		}
		if (i*i == n) if (m%i == 0) res.push_back(i);
		sort(res.begin(), res.end());
		return res;
	}




}



