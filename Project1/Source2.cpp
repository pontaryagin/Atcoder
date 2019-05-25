#include "MyHeader.h"
//#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"





// ============================ Header  =================================



int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll q;
	cin >> q;
	ll sumb = 0;
	ll sumal = 0;
	ll sumar = 0;

	set<ll> seta;
	priority_queue<ll> pql;
	pq_greater<ll> pqr;

	ll cent = 0;
	rep(i, 0, q) {
		ll flag;
		cin >> flag;
		if (flag == 1) {
			ll a, b;
			cin >> a >> b;
			sumb += b;
			if (i == 0) {
				pql.push(a);
				cent = a;
				sumal += a;

			}
			else if (pqr.size() == pql.size()) {
				ll l = pql.top();
				ll r = pqr.top();
				if (a <=l) {
					cent = max(a, l);
					pql.push(a);
					sumal += a;
				}
				else{
					pqr.push(a);
					sumar += a;
					cent = pqr.top();
					pqr.pop(); pql.push(cent);
					sumar -= cent; sumal += cent;

				}
			}
			else {
				ll l = pql.top();
				ll r = (pqr.size()==0 ?0:pqr.top());
				if (a < cent) {
					pql.push(a);
					sumal += a; sumal -= pql.top(); sumar += pql.top();
					pqr.push(pql.top()); pql.pop();
					cent = pql.top();

				}
				else {
					pqr.push(a);
					sumar += a;

				}
			}


		}
		else {
			cout << cent << " " << sumb + SZ(pql) * cent - sumal + sumar - SZ(pqr)*cent <<endl;
		}
	}
	
	return 0;

}
