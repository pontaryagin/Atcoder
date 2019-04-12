#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
//#include "NumberTheory.h"
//#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"

// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);

	ll x, y, z, k;
	cin >> x >> y >> z >> k;
	vll a(x), b(y), c(z),a2(x),b2(y),c2(y);
	rep(i, 0, x) {
		cin >> a[i];
	}
	rep(i, 0, y) {
		cin >> b[i];
	}
	rep(i, 0, z) {
		cin >> c[i];
	}

	sort(all(a)), sort(all(b)), sort(all(c));
	rep(i, 0, x) {
		a2[i] = i == 0 ? 0 : a[i] - a[i - 1];
	}
	rep(i, 0, y) {
		b2[i] = i == 0 ? 0 : b[i] - b[i - 1];
	}
	rep(i, 0, z) {
		c2[i] = i == 0 ? 0 : c[i] - c[i - 1];
	}
	ll cnt = 1;
	ll sum = a[0]+b[0]+c[0];
	ll i_a=1, i_b=1,i_c=1;
	while (cnt < k) {
		ll m = min(a2[i_a], min(b2[i_b], c2[i_c]));
		if (a2[i_a] == m) {
			i_a++;
			rep (i, 0, i_b) {
				sum += b
			}

		}
	}



	return 0;
	
}


