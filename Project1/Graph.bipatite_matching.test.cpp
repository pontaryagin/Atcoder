#define PROBLEM http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_7_A&lang=jp

#include "MyHeader.h"
#include "Graph.h"

int main() {

	ll x, y, m; cin >> x >> y >> m;
	BipartiteMatching g(x, y);
	rep(i, 0, m) {
		ll x, y; cin >> x >> y;
		g.push(x, y);
	}
	cout << g.bipartite_match() << endl;
	return 0;
}