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

			}
			else if (pqr.size() == pql.size()) {
				ll l = pql.top();
				ll r = pqr.top();
				if (a <=l) {
					cent = max(a, l);
					
					suma -= pqr.size() * abs(cent - l);
					pql.push(a);
				}
				else{
					cent = min(r,a);	
					suma += abs(a - cent);
					suma += (pql.size() - 1) * abs(cent - l);
					suma -= (pqr.size() + 1) * abs(cent - l);
					pqr.pop(); pql.push(cent); pqr.push(a);

				}
			}
			else {
				ll l = pql.top();
				ll r = (pqr.size()==0 ?0:pqr.top());
				if (a < cent) {
					pql.push(a);
					ll l = pql.top(); pql.pop();
					cent = pql.top();
					pqr.push(l);
					suma += abs(a - l);
					suma -= pql.size() * abs(cent - l);
					suma += pqr.size() * abs(cent - l);

				}
				else {
					pqr.push(a);
					suma += abs(a - cent);

				}
			}


		}
		else {
			cout << cent << " " << sumb + suma <<endl;
		}
	}
	
	return 0;

}
