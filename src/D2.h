#include "MyHeader.h"
#include <optional>

struct D2 {
    enum Dir {U, D, L, R};
    static const inline vector<Dir> Dirs = { U, D, L, R };
    D2(ll h, ll w) : h(h), w(w) {};
    bool in(ll n, Dir d, ll k = 1) {
        auto [H, W] = (*this)(n);
        switch (d)
        {
        case U: return (H + k < h) && (H + k >=0);
        case D: return (H - k < h) && (H - k >= 0);
        case R: return (W + k >= 0) && (W + k < w );
        case L: return (W - k >= 0) && (W - k < w);
        default: throw invalid_argument("not direction");
        }
    };
    optional<ll> next(ll n, Dir d, ll k = 1) {
        if (!in(n, d, k)) return nullopt;
        switch (d)
        {
        case U: return n + w * k;
        case D: return n + -w * k;
        case L: return n + -k;
        case R: return n + k;
        default: throw invalid_argument("not direction");
        }
    }
    pair<ll, ll> operator()(ll val) {
        return { val / w, val % w };
    }
    ll operator()(ll h_, ll w_) {
        return h_ * w + w_;
    }
    ll h, w;
};
