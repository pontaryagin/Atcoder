#include "MyHeader.h"
#include "Graph.h"
#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"


// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll n, k; cin >> n >> k;
	Graph g(n);
	vll a(n - 1), b(n - 1);
	rep(i, 0, n - 1) {
		cin >> a[i] >> b[i];
		a[i]--; b[i]--;
		g.push({ a[i],b[i] });
		g.push({ b[i],a[i] });
	}
	Combination cmb;
	GraphDFS dfs(g);
	modint<> res = k;
	ll sz =  g.out(0).size(); 
	if (sz > k - 1) {
		res = 0;
	}
	else {
		res *= cmb.P(k-1, sz);
	}
	dfs([&](Edge & e)
		{
			
			if (g.out(e.to).size() <= k - 1) {
				auto sz = g.out(e.to).size();
				if (sz == 1) {
					return;
				}
				sz -= 1;
				ll k_ = k - 2;
				res *= cmb.P(k_, sz);

			}
			else
				res = 0;
		}
	, 0);

	cout << res << endl;

	return 0;

}