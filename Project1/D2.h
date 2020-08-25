#include "MyHeader.h"

struct D2 {
	enum Dir {U, D, L, R};
	static inline vector<Dir> Dirs = { U, D, L, R };
	D2(ll h, ll w) : h(h), w(w) {};
	bool in(ll n, Dir d, ll k = 1) {
		switch (d)
		{
		case U: return (H(n) + k < h) && (H(n) + k >=0);
		case D: return (H(n) - k < h) && (H(n) - k >= 0);
		case R: return (W(n) + k >= 0) && (W(n) + k < w );
		case L: return (W(n) - k >= 0) && (W(n) - k < w);
		default: throw invalid_argument("not direction");
		}
	};
	ll next(ll n, Dir d, ll k = 1) {
		switch (d)
		{
		case U: return n + w * k;
		case D: return n + -w * k;
		case L: return n + -k;
		case R: return n + k;
		default: throw invalid_argument("not direction");
		}
	}
	ll H(ll n) { return n / w; }
	ll W(ll n) { return n % w; }
	ll h, w;
};


