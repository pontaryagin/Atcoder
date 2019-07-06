#include "MyHeader.h"
//#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"


// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);
	
	//E - Range Minimum Queries
	ll n, k , q; cin >> n>>k>>q;
	vll a(n);
	read(a);
	
	LazySegmentTree<M::min_indexed_t<>> seg(n);
	rep(i, 0, n) seg.update(i, i + 1, { a[i],i });

	ll res = 
	rep(i, 0, n - q + 1)
	{

	}


	return 0;

}