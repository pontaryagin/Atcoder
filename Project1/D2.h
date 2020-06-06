#include "MyHeader.h"

struct D2 {
	D2(ll h, ll w) : h(h), w(w), Dir{1, -1, w, -w} {};
	bool in(ll n, ll d) {
		if (n % w == 0 && d == -1)
			return false;
		if (n % w == w - 1 && d == 1)
			return false;
		if (n / w == 0 && d == -w)
			return false;
		if (n / w == h - 1 && d == w)
			return false;
		return true;
	};
	ll h, w;
	vll Dir;
};