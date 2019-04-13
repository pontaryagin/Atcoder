#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
//#include "NumberTheory.h"
//#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"

// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);

	queue<ll> q;
	ll n, k;
	cin >> n >> k;
	string S;
	cin >> S;
	ll m = 0;
	q.push(0);
	ll zero = 0;
	if (S[0] == '0') {
		zero++;
	}
	rep(i, 1, n) {
		if (zero == k) {
			if (S[i-1] - '0' == 1 && S[i ] - '0' == 0) {

				ll x = q.front(); q.pop();
				m = max(m, i - x);

			}
			else if (S[i - 1] - '0' == 0 && S[i] - '0' ==1) {
				q.push(i);


			}
		}
		else {
			if (S[i - 1] - '0' == 1 && S[i] - '0' == 0) {
				k--;
			}
			else if (S[i - 1] - '0' == 0 && S[i] - '0' == 1) {
				q.push(i);
			}
		}

	}
	m = max(m, n - q.front());
	cout << m;

	return 0;
	
}

//
//ll T;
//cin >> T;
//vll res(T);
//rep(k, 0, T) {
//	ll n;
//	cin >> n;
//	vll W(n);
//	rep(i, 0, n) {
//		cin >> W[i];
//	}
//
//}
//
//rep(k, 0, T) {
//	cout << "Case #" << k << ": " << res[k] << endl;
//}