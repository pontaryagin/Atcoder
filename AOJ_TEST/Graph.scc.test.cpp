#define PROBLEM "http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_3_C&lang=jp"
#include "Graph.h"

int main() {
    ll n, m; cin >> n >> m;
    Digraph g(n);
    rep(i, 0, m) {
        ll s, t; cin >> s >> t;
        g.push({ s,t });
    }
    vll par(n);
    Digraph dag = g.scc(par);

    ll q; cin >> q;
    rep(i, 0, q) {
        ll u, v; cin >> u >> v;
        cout << (par[u] == par[v] ? 1 : 0) << endl;
    }
    return 0;
}
