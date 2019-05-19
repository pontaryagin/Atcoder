#include "MyHeader.h"
//#include "Graph.h"
#include "NumberTheory.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"





// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	vec_t<1,pair<ll, pll>> x(2);
	read_v(x);


	ll m, k;
	cin >> m >> k;
	if (k>(1<<m)-1 ) {
		cout << -1 << endl;
	}
	else if (k != 0) {
		if (m == 1) {
			if (k == 0) cout<< "0 0 1 1" <<endl;
			else cout<<-1 <<endl;
			return 0;
		}
		vll res;
		res.push_back(0); res.push_back(k); res.push_back(0);
		rep(i, 1, 1 << m) {
			if(i!= k)
				res.push_back(i);
		}
		res.push_back(k);
		rrep(i, 1, 1 << m) {
			if (i != k)
				res.push_back(i);
		}
		write_v(res);
	}
	else {
		vll res = seq(0, 1 << m);
		rep(i, 0, res.size()) {
			if (i != res.size() - 1)
				cout << res[i] << " " << res[i] << " ";
			else
				cout << res[i] << " " << res[i] << endl;

		}
	}

	return 0;

}