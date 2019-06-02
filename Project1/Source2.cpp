#include "MyHeader.h"
//#include "Graph.h"
#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"





// ============================ Header  =================================

template<typename Inputs>
void sort_uniq(Inputs& inputs) {
	sort(all(inputs));
	inputs.erase(unique(all(inputs)), inputs.end());
}

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll q;
	cin >> q;
	vector<modint<>> x(q), d(q);vll n(q);
	vector<modint<>> S(MOD+1, 1);
	rep(i, 1, MOD) {
		S[i + 1] = S[i] * modint<>(i);
	}
	rep(i, 0, q) {
		cin >> x[i] >> d[i] >> n[i];
		if (n[i] > MOD)cout << 0 << endl;
		else {
			if (d[i] == 0) {
				cout << POW(x[i] ,n[i] % MOD)<<endl;
			}
			else if (n[i] == 1) {
				cout << x[i].a << endl;
			}
			else {
				auto D = x[i]/d[i];
				modint<> res =1 ;
				if (D.a + n[i] - 1 < MOD) {
					res *= S[D.a + n[i] ] /S[D.a];

				}
				else {
					res = 0;
				}
				res *= POW(d[i], n[i]%MOD - 2);
				cout << res.a << endl;
			}
		}

	}





	





	return 0;

}
