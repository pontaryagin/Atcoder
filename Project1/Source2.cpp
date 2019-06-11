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
	

	//Int l, a, b;
	//cin >> l >> a >> b >> MOD;

	//using T = tuple<Int, Int, Int>;
	//vector<T> vt;
	//Int c = a + b * (l - 1);

	//auto calcLen = [](Int x) {return (Int)to_string(x).size(); };
	//Int x = a;
	//while (x <= c) {
	//	Int L = (x - a) / b;
	//	Int R = l;
	//	if (calcLen(x) == calcLen(c)) {
	//		vt.emplace_back(x, L, R);
	//		break;
	//	}
	//	while (L + 1 < R) {
	//		Int M = (L + R) >> 1;
	//		if (calcLen(a + b * M) == calcLen(x)) L = M;
	//		else R = M;
	//	}
	//	vt.emplace_back(x, (x - a) / b, R);
	//	x = a + b * R;
	//}
	//reverse(vt.begin(), vt.end());
	//using M = Mint<Int>;
	//using SM = SquareMatrix<2, M>;
	//M ans{ 0 };
	//M mul{ 1 };
	//for (auto st : vt) {
	//	Int x, L, R;
	//	tie(x, L, R) = st;
	//	Int len = calcLen(x);

	//	auto makeB =
	//		[&](Int num) {
	//		SM A;
	//		A[0][0] = M(1); A[0][1] = M(1);
	//		A[1][0] = M(0); A[1][1] = M(10).pow(len);
	//		return M(b)* (A.pow(num))[0][1];
	//	};

	//	using P = pair<M, Int>;
	//	auto merge =
	//		[&](P a, P b) {
	//		M x = a.first;
	//		M y = b.first;
	//		Int p = a.second;
	//		Int q = b.second;

	//		M z{ 0 };
	//		z += x * (M(10).pow(len)).pow(q);
	//		z += y + M(p) * makeB(q);

	//		return P(z, p + q);
	//	};

	//	Int n = R - L;
	//	P add(0, 0);
	//	P res(x, 1);
	//	while (n) {
	//		if (n & 1) add = merge(add, res);
	//		res = merge(res, res);
	//		n >>= 1;
	//	}
	//	//cout<<add.first.v<<" "<<add.second<<endl;
	//	ans += add.first * mul;
	//	mul *= (M(10).pow(len)).pow(R - L);
	//}

	//cout << ans.v << endl;


	return 0;

}