
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
	
	string s, t;
	cin >> s >> t;
	auto dp = make_v(s.size() + 1, t.size() + 1);
	//auto dps = make_v<string>(s.size() + 1, t.size() + 1);
	rep(i, 0, s.size())rep(j, 0, t.size()) {
		dp[i + 1][j + 1] = max(dp[i + 1][j], dp[i][j + 1]);
		//dps[i + 1][j + 1] = (dp[i + 1][j] > dp[i][j + 1] ? dps[i + 1][j] : dps[i][j + 1]);
		if (s[i] == t[j]) {
			if (dp[i][j] + 1 > max(dp[i + 1][j], dp[i][j + 1])){
				dp[i + 1][j + 1] = dp[i][j] + 1;
				//dps[i + 1][j + 1] = dps[i][j] + s[i];
			}
		}
		//dps[i][j].clear();
		//dps[i][j].shrink_to_fit();
	}
	string res;
	ll i = s.size(); ll j = t.size();
	while (dp[i][j] > 0) {
		if (dp[i][j] == dp[i - 1][j]) {
			--i;
		}
		else if (dp[i][j - 1] == dp[i][j]) {
			--j;
		}
		else {
			res.push_back(s[i-1]);
			--i; --j;
		}
	}
	reverse(all(res));
	cout << res << endl;
	return 0;

}
