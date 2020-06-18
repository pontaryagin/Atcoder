#include "MyHeader.h"

struct D2 {
	enum Dir {U, D, L, R};
	static inline const Dir Dirs[] = { U, D, L, R };
	D2(ll h, ll w) : h(h), w(w) {};
	bool in(ll n, D2::Dir d, ll k = 1) {
		if ((W(n) <= k-1) && d == L)
			return false;
		if ((W(n) >= w - 1 - (k-1)) && d == R)
			return false;
		if ((H(n) <= k-1) && d == D)
			return false;
		if ((H(n) >= h - 1 - (k - 1)) && d == U)
			return false;
		return true;
	};
	ll next(Dir d, ll k = 1) {
		switch (d)
		{
		case U: return w * k;
		case D: return -w * k;
		case L: return -k;
		case R: return k;
		default: throw invalid_argument("not direction");
		}
	}
	ll H(ll n) { return n / w; }
	ll W(ll n) { return n % w; }

	ll h, w;
};


