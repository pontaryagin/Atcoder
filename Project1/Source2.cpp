#include "MyHeader.h"
//#include "Graph.h"
#include "NumberTheory.h"
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

	ll n;
	modint<> A, B, C;
	cin >> n >> A>> B >> C;
	auto a = A/(A+B) ; auto b = B/(A+B) ;

	modint<> res = 0;
	auto cmb = Combination( MOD, 200000 );
	rep(m, n, 2 * n) {
		//cout << cmb(m - 1, n - 1) << endl;;
		//cout << a * a << endl;
		//cout << POW(b, m - n) << endl;;
		res += modint<>(cmb(m - 1, n - 1)) * POW(a, n) * POW(b, m - n)*m;
		//cout <<"res"<< res << endl;
		res += modint<>(cmb(m - 1, n - 1)) * POW(b, n) * POW(a, m - n)*m;
		//cout << "res" << res << endl;
	}
	cout << (res*(100)/(A+B)) << endl;
	return 0;

}
