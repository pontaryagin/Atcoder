#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
#include "NumberTheory.h"
#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"



// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll N;
	cin >> N;
	vector<double> p(N);
	rep(i, 0, N)cin >> p[i];
	auto dp = make_v<double>(N + 1, N + 1);
	dp[0][0] = 1;
	rep(i, 0, N)rep(j, 0, N ) {
		dp[i + 1][j + 1] += dp[i][j]*p[j];
		dp[i][j+1] += dp[i][j] * (1-p[j]);
	}
	double res = 0;
	rep(i, 0, N + 1) {
		if (i > N / 2)
			res += dp[i][N];
	}
	cout << res << endl;
	return 0;

}
