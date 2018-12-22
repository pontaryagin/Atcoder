
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
// #include "bits/stdc++.h"

using namespace std;
//debug

#define rep(i, N) for (int i = 0; i < N; i++)

#define pb push_back

typedef long long ll;
typedef pair<ll, ll> pll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<string> vs;



#define all(a)  (a).begin(),(a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define vec(a) vector<a>


int _Dijk() {
	ll n;
	cin >> n;
	vector<vector<pll > > adj(n);//adj:隣接リスト。adj[NodeFrom][j]=(wight,NodeTo)
	ll u, k, v, c;
	rep(i, n) {
		cin >> u >> k;
		adj[i].resize(k);
		rep(j, k) {
			cin >> v >> c;
			adj[i][j] = make_pair(c, v);
		}
	}
	ll s = 0;
	//dijkstra

	vector<ll> res(n, 1 << 30);//store result
	vector<int> state(n, 0);//store state. 0=>imcomplete. 1=>complete.

	priority_queue<pll, vector<pll>, greater<pll> > pq;//priority_queueはデフォルトで昇順なので降順に直す
	pq.push({ 0,s });
	while (!pq.empty()) {

		pll top = pq.top(); pq.pop();//get min and pop it.
		ll newNode = top.second;

		if (state[newNode] != 0) continue; //if newNode has already been visited then skip.
		res[newNode] = top.first;
		state[newNode] = 1;

		for (pll pr : adj[newNode]) {
			if (state[pr.second] == 0) {
				pq.push({ top.first+pr.first,pr.second });
			}
		}
	}

	rep(i, n) {
		cout << i << " " << res[i] << endl;
	}
	return 0;
}

