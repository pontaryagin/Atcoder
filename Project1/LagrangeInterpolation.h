#include "NumberTheory.h"

//using modt = modint<-11232>;
//modt::set_modulo(p);
//vector<modt> a(p); read(a);
//
//vector<vector<modt>> l(p, vector<modt>(p + 1)); // Lagrange Polynomial
//Combination cmb(p);
//vector<vector<modt>> m(p + 1, vector<modt>(p + 1));
//// m[i] = x * (x-1) * (x-2) * ... * (x-i);
//m[0][0] = 1;
//rep(i, 0, p) {
//	rep(j, 0, p + 1) {
//		m[i + 1][j] = (j > 0 ? m[i][j - 1] : 0) + m[i][j] * (-i);
//	}
//}
//auto& M = m[p];
//rep(i, 0, p) {
//	rrep(j, 0, p) {
//		l[i][j] = l[i][j + 1] * i + M[j + 1];
//	}
//}
//
//rep(i, 0, p) {
//	rrep(j, 0, p) {
//		l[i][j] *= cmb.FacInv(i);
//		l[i][j] *= cmb.FacInv(p - 1 - i);
//		if ((p - 1 - i) % 2 == 1) l[i][j] *= -1;
//	}
//}
//vector<modt> res(p);
//rep(i, 0, p) {
//	rrep(j, 0, p) {
//		res[j] += a[i] * l[i][j];
//	}
//}