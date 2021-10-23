#define PROBLEM "http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DPL_1_D"

#include "MyHeader.h"
//#include "NumberTheory.h"
#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "FFT.h"
#include "Graph.h"
//#include "LagrangeInterpolation.h"
//#include "UnionFind.h"
//#include "Text.h"
//#include "Convolution.h"
//#include "Timer.h"


int main() {

    cin.tie(0);
    cout.tie(0);
    ios::sync_with_stdio(false);
    cout << fixed << setprecision(12);


    ll n; cin>> n;
    vll a(n); cin >> a;
    vll a2 = a; sort_uniq(a2);
    auto a_ind = inv_map(a2);
    rep(i, 0, n) {
        a[i] = a_ind[a[i]];
    }
    segment_tree<M::max_t<ll>> seg(n, 0LL);
    ll res = 0;
    rep(i, 0, n) {
        ll v = seg.query(0, a[i]);
        seg.update(a[i], v + 1);
        chmax(res, v + 1);

    }

    cout << res<< endl;
    return 0;
}
