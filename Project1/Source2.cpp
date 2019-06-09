#include "MyHeader.h"
#include "Graph.h"
//#include "NumberTheory.h"
#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"





// ============================ Header  =================================



int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll n, m; cin >> n >> m;
	vpll ab(m);
	read(ab);

	UnionFind uf(n);
	Graph g(n);
	rep(i, 0, m) {
		uf.unite(ab[i].fir-1, ab[i].sec-1);
		g.push(Edge{ ab[i].first - 1,ab[i].sec - 1 });

	}
	bool ok = 1;
	rep(i, 0, n) {
		if ((g.out(i).size() + g.in(i).size()) % 2) {
			ok = 0;
		}
	}
	if (ok)
		cout << "YES" << endl;
	else
		cout << "NO" << endl;

	return 0;

}