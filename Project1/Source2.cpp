#include "MyHeader.h"
#include "Graph.h"
//#include "NumberTheory.h"
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
	vll x(q), d(q), n(q);
	read(x, d, n);

	





	return 0;

}
