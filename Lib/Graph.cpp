
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
#include <time.h>
//#include "boost/variant.hpp"

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

template<typename T1, typename T2> void chmin(T1 &a, T2 b) { if (a > b) a = b; }
template<typename T1, typename T2> void chmax(T1 &a, T2 b) { if (a < b) a = b; }




vll get_topologically_sorted_nodes(const vvll& graph)
{
	// graph needs to be represented by adjacent list.
	ll nodeSize = graph.size();

	// find root
	vll roots;
	vll inDegree(nodeSize);
	rep(i, 0, nodeSize)
	{
		for (ll sibling : graph[i]) {
			inDegree[sibling]++;
		}
	}

	rep(i, 0, nodeSize) {
		if (inDegree[i] == 0) {
			roots.push_back(i);
		}
	}

	stack<ll> parents;
	for (ll i : roots)
		parents.push(i);

	vll sortedNodes;
	while (!parents.empty()) {
		ll parent = parents.top();
		parents.pop();
		sortedNodes.push_back(parent);
		for (ll sibling : graph[parent]) {
			inDegree[sibling]--;
			if (inDegree[sibling] == 0) {
				parents.push(sibling);
			}
		}
	}
	return sortedNodes;
}


