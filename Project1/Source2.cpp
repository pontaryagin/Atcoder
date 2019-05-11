#include "MyHeader.h"
#include "Graph.h"
//#include "NumberTheory.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"






// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);


	ll n;
	cin >> n;
	pq_greater<Edge> es;
	vll p(n - 1);
	vll h(n - 1);
	rep(i, 0, n - 1) {
		Edge e;

		cin >> e.from;
		p[i] = e.from;
		cin >> e.cost;
		h[i] = e.cost;
		e.to = i+1;
		es.emplace(move(e));
	}
	vll removed(n);
	ll res = 0;
	function<void(ll,ll,ll)> rec = [&](ll i, ll cost, ll res_ ) {
		if (removed[i])
			return;
		ll j = p[i-1];
		ll c = h[i-1];
		c -= cost;
		if (c == 0) {
			res_++;
			removed[i] = 1;
		}
		else
			es.emplace(Edge(j, i, c));

		if (j != 0) {
			rec(j, cost, res_);
		}
		else {
			res += res_;
		}
	};
	while(es.size()) {

		
		auto e = es.top(); es.pop();
		rec(e.to, e.cost, 0);
		
	}

	cout << res+1 << endl;



	return 0;

}

