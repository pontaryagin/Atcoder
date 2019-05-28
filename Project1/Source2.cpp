#include "MyHeader.h"
//#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
#include "SegmentTree.h"
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

	
	int n, q;
	cin >> n >> q;
	LazySegmentTree<Monoid::sum_t<>> ch(n);
	for (int i = 0; i < q; i++) {
		int c, s, t, x;
		cin >> c;
		if (c) {
			cin >> s >> t;
			cout << ch.query(s - 1, t) << endl;
		}
		else {
			cin >> s >> t >> x;
			ch.update(s - 1, t, x);
		}
	}

	
	return 0;

}
