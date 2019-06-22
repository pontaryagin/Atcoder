#include "MyHeader.h"
//#include "Graph.h"
#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"
#include "LinearAlgebra.h"

// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);
	

	vll len(37,INF);
	ll l, a, b, m;
	cin >> l >> a >> b >> m;
	ll cursz = 0;
	rep(i, 0, l) {
		ll sz = to_string(a + b * i).size();
		if (sz != cursz) {
			len[sz = i];
			cursz = sz;
		}
	
	}

	rep(r, 0, SZ(len)) {
		ll lenth = min(len[r + 1], l) - len[r];
		SquareMatrix<3, modint<>> sq;
		sq.dat[0] = { POW(modint<>(10),lenth),1,1 };
		sq.dat[1] = { 0,1,1 };
		sq.dat[2] = { 0,0,1 };
		sq = sq.pow(5);
		sq =sq* sq;


	}



	return 0;

}