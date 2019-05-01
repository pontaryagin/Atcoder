
#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
#include "NumberTheory.h"
//#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"



// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	
	ll N, W;
	cin >> N>>W;
	vll w(N), v(N);
	rep(i, 0, N)cin >> w[i] >> v[i];
	auto dp = make_v(N+1, W + 1);
	fill_v(dp, -1);
	dp[0][0] = 0;
	rep(i, 0, N){
		dp[i + 1] = dp[i];
		rep(wei, 0, W + 1) {
			if (wei + w[i] <= W && dp[i][wei] != -1) {
				chmax(dp[i + 1][wei + w[i]], dp[i][wei] + v[i]);
			}
		}
	}
	cout << *max_element(all(dp[N]))<<endl;
	return 0;

}
