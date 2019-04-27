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

	ll T;
	cin >> T;
	rep(test, 1, T + 1) {
		[&] {
			ll r, c, h, v;
			cin >> r >> c >> h >> v;
			auto choco = make_v<ll>(r, c);
			auto S = make_v<ll>(r + 1, c + 1);

			rep(i, 0, r) {
				//vll ss(c + 1);
				rep(j, 0, c) {
					char c;
					cin >> c;
					choco[i][j] = c == '@';
					//ss[j + 1] = ss[j] + int(choco[i][j] == '@');
					//S[i + 1][j + 1] = S[i][j + 1] + ss[j + 1];
				}
			}
			S = commulativeSum(choco);
			vll rows;
			vll columns;

			ll divn = (h + 1) * (v + 1);
			if (S[r][c] % (divn) != 0)
			{
				cout << "Case #" << test << ": IMPOSSIBLE" << endl;
				return;
			}
			ll unit = S[r][c] / divn;
			rep(i, 1, r) {
				if ((rows.size() + 1) * (v + 1) * unit == S[i][c]) {
					rows.push_back(i);
					if (rows.size() == h)
						break;
				}
			}
			rep(i, 1, c) {
				if ((columns.size() + 1) * (h + 1) * unit == S[r][i]) {
					columns.push_back(i);
					if (columns.size() == v)
						break;
				}
			}
			if (rows.size() != h || columns.size() != v) {
				cout << "Case #" << test << ": IMPOSSIBLE" << endl;
				return;
			}
			// check
			rep(i, 0, h)rep(j, 0, v) {
				if (S[rows[i]][columns[j]] != unit * (j + 1) * (i + 1)) {
					cout << "Case #" << test << ": IMPOSSIBLE" << endl;
					return;

				}


			}
			cout << "Case #" << test << ": POSSIBLE" << endl;
		}();
	}

	return 0;

}
