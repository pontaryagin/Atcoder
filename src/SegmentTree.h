#pragma once
#include "MyHeader.h"


template<typename T, typename S = nullptr_t>
struct has_unit
    :public false_type
{
};

template<typename T>
struct has_unit<T, typename conditional<false, decltype(T::unit()), nullptr_t>::type >
    :public true_type
{
};
template<typename T, typename S = nullptr_t>
struct has_append
    :public false_type
{
};

template<typename T>
struct has_append<T, typename conditional<false, decltype(T::append(T::underlying_type(), T::underlying_type())), nullptr_t>::type >
    :public true_type
{
};



template <typename Monoid>
struct segment_tree
{

    using underlying_type = typename  Monoid::underlying_type;

    segment_tree(ll a_n, underlying_type unit = Monoid::unit()) : size_original(a_n), unit(unit)
    {
        vector<underlying_type> initial_value = vector<underlying_type>(a_n, unit);
        segment_tree_impl(a_n, initial_value);
    }

    segment_tree(ll a_n, vector<underlying_type>& initial_value, underlying_type unit = Monoid::unit()) : size_original(a_n), unit(unit)
    {
        segment_tree_impl(a_n, initial_value);
    }

    void update(int i, underlying_type z) { // 0-based
        assert(0 <= i && i < 2 * n - 1);
        a[i + n - 1] = z;
        for (i = (i + n) / 2; i > 0; i /= 2) { // 1-based
            a[i - 1] = Monoid::append(a[2 * i - 1], a[2 * i]);
        }
    }

    underlying_type query(ll l, ll r) { // 0-based, [l, r)
        underlying_type lacc = unit, racc = unit;
        assert(l <= r && r <= n);
        l += n; r += n;
        for (; l < r; l /= 2, r /= 2) { // 1-based loop, 2x faster than recursion
            if (l % 2 == 1) lacc = Monoid::append(lacc, a[(l++) - 1]);
            if (r % 2 == 1) racc = Monoid::append(a[(--r) - 1], racc);
        }
        return Monoid::append(lacc, racc);
    }

    underlying_type query(ll i) { // return value at i
        assert(0 <= i && i < size_original);
        return a[i + n - 1];
    }

    ll size() { return size_original; }

private:
    ll size_original;
    ll n;
    vector<underlying_type> a;
    underlying_type unit;
    void segment_tree_impl(ll a_n, vector<underlying_type>& initial_value)
    {
        assert(a_n == initial_value.size());
        n = 1; while (n < a_n) n *= 2;
        a.resize(2 * n - 1, unit);
        rep(i, 0, initial_value.size()) {
            a[i + (n - 1)] = initial_value[i];
        }
        rrep(i, 0, n - 1) a[i] = Monoid::append(a[2 * i + 1], a[2 * i + 2]); // propagate initial values
    }


};


// template <typename Monoid>
using Monoid = M::max_t<ll>;
struct segment_tree_2d
{

    using underlying_type = typename  Monoid::underlying_type;
    using data_type = vector<vector<underlying_type>>;

    segment_tree_2d(ll h, ll w, underlying_type unit = Monoid::unit()) 
        : unit(unit)
    {
        auto initial_value = data_type(h, vector(w, unit));
        segment_tree_impl(move(initial_value));
    }

    segment_tree_2d(data_type initial_value, underlying_type unit = Monoid::unit())
        : unit(unit)
    {
        segment_tree_impl(move(initial_value));
    }

    void update(int hh, int ww, underlying_type z) { // 0-based
        assert(0 <= hh && hh < h_org && 0<= ww && ww < w_org);
        a[_base_pos1(hh)][_base_pos2(ww)] = z;
        // rep (i = (i + n) / 2; i > 0; i /= 2) { // 1-based
        auto base_h = _base_pos1(hh);
        auto base_w = _base_pos2(ww);
        // {base_h} * [0, 2*w-1) is updated
        for (auto _w = _par(base_w); _w >= 0; _w = _par(_w)){
            a[base_h][_w] = Monoid::append(a[base_h][_l(_w)], a[base_h][_r(_w)]);            
        }
        for (auto _w = base_w; _w >= 0; _w = _par(_w)){
            for (auto _h = _par(base_h); _h >= 0; _h = _par(_h)){
                a[_h][_w] = Monoid::append(a[_l(_h)][_w], a[_r(_h)][_w]);            
            }
        }
    }

    underlying_type query(ll hh1, ll ww1, ll hh2, ll ww2) { // 0-based, [l, r)
        underlying_type lacc = unit, racc = unit;
        assert(l <= r && r <= n);
        l += n; r += n;
        for (; l < r; l /= 2, r /= 2) { // 1-based loop, 2x faster than recursion
            if (l % 2 == 1) lacc = Monoid::append(lacc, a[(l++) - 1]);
            if (r % 2 == 1) racc = Monoid::append(a[(--r) - 1], racc);
        }
        return Monoid::append(lacc, racc);
    }

    underlying_type query_h(ll hh, ll ww1, ll ww2) { // 0-based, [l, r)
        underlying_type lacc = unit, racc = unit;
        assert(l <= r && r <= n);
        l += n; r += n;
        for (; l < r; l /= 2, r /= 2) { // 1-based loop, 2x faster than recursion
            if (l % 2 == 1) lacc = Monoid::append(lacc, a[(l++) - 1]);
            if (r % 2 == 1) racc = Monoid::append(a[(--r) - 1], racc);
        }
        return Monoid::append(lacc, racc);
    }

    underlying_type query(ll hh, ll ww) { // return value at i
        assert(0 <= hh && hh < h_org && 0<= ww && ww < w_org);
        return a[_base_pos1(hh)][_base_pos2(ww)];
    }

    pll size() { return {h_org, w_org}; }

private:
    ll h_org, w_org;
    ll h, w;
    data_type a;
    underlying_type unit;
    int _base_pos1(int i) { return i + h-1; }
    int _base_pos2(int i) { return i + w-1; }
    int _l(int i) { return 2*i+1; }
    int _r(int i) { return 2*i+2; }
    int _par(int i) { return (i-1)/2; }

    void segment_tree_impl(data_type initial_value)
    {
        h_org = initial_value.size();
        w_org = initial_value.at(0).size();
        h = 1; while (h < h_org) h *= 2;
        w = 1; while (w < w_org) w *= 2;
        a = move(initial_value);
        a.resize(2 * h - 1, vector<underlying_type>(2*w-1, unit));
        // [h-1, 2*h-1) * [w-1, 2*w-1) is initialized here
        rep(i, 0, initial_value.size()) {
            rep(j,0,initial_value[i].size()){
                a[_base_pos1(i)][_base_pos2(j)] = initial_value[i][j];
            }
        }
        // propagate w.r.t ww for fixed hh
        // [h-1, 2*h-1) * [0, 2*w-1) is initialized here
        rrep(hh, 0, h) {
            rrep(ww,0, w-1){
                auto base_h = _base_pos1(hh);
                a[base_h][ww] = Monoid::append(a[base_h][_l(ww)], a[base_h][_r(ww)]);
            }
        }
        // propagate w.r.t hh for fixed ww
        // [0, 2*h-1) * [0, 2*w-1) is initialized here
        rrep(ww,0, 2*w-1){
            rrep(hh, 0, h-1) {
                a[hh][ww] = Monoid::append(a[_l(hh)][ww], a[_r(hh)][ww]);
            }
        }
    }


};



namespace M {

    template <typename T = ll>
    struct sum_t {
        typedef T underlying_type;
        static underlying_type unit() { return 0; }
        static underlying_type append(underlying_type a, underlying_type b) { return a + b; }
        static underlying_type iterate(underlying_type a, int n) { return a * n; }
    };
    template <typename T = ll>
    struct xor_t {
        typedef T underlying_type;
        static underlying_type unit() { return 0; }
        static underlying_type append(underlying_type a, underlying_type b) { return a ^ b; }
    };
    template <typename T = ll>
    struct or_t {
        typedef T underlying_type;
        static underlying_type unit() { return 0; }
        static underlying_type append(underlying_type a, underlying_type b) { return a | b; }
    };
    template <typename T = ll>
    struct gcd_t {
        typedef T underlying_type;
        static underlying_type unit() { return INF; }
        static underlying_type append(underlying_type a, underlying_type b) { return gcd(a, b); }
    };

    template<typename S, typename T>
    struct pair_t {
        typedef pair<typename S::underlying_type, typename T::underlying_type> underlying_type;
        static underlying_type unit() { return make_pair(S::unit(), T::unit()); }
        static underlying_type append(underlying_type a, underlying_type b) { return make_pair(S::append(a.first, b.first), T::append(a.second, b.second)); }
        static underlying_type iterate(underlying_type a, int n) { return make_pair(S::iterate(a.first, n), T::iterate(a.second, n)); }
    };

    template <typename T = ll>
    struct min_t {
        typedef T underlying_type;
        static underlying_type unit() { return numeric_limits<T>::max(); }
        static underlying_type append(underlying_type a, underlying_type b) { return min(a, b); }
        static underlying_type iterate(underlying_type a, size_t n) { return a; }
    };

    template <typename T = ll>
    struct max_t {
        typedef T underlying_type;
        static underlying_type unit() { return numeric_limits<T>::min(); }
        static underlying_type append(underlying_type a, underlying_type b) { return max(a, b); }
        static underlying_type iterate(underlying_type a, size_t n) { return a; }
    };

    template <typename T = ll, typename IndexType = ll>
    struct min_indexed_t {
        typedef pair<typename min_t<T>::underlying_type, IndexType> underlying_type;
        static underlying_type unit() { return make_pair(numeric_limits<T>::max(), IndexType{}); }
        static underlying_type append(underlying_type a, underlying_type b) { return min(a, b); }
        static underlying_type iterate(underlying_type a, int n) { return a; }
    };
    template <typename T = ll, typename IndexType = ll>
    struct max_indexed_t {
        typedef pair<typename min_t<T>::underlying_type, IndexType> underlying_type;
        static underlying_type unit() { return make_pair(numeric_limits<T>::min(), IndexType{}); }
        static underlying_type append(underlying_type a, underlying_type b) { return max(a, b); }
        static underlying_type iterate(underlying_type a, int n) { return a; }
    };

    struct linear_t {
        typedef pd underlying_type;
        static underlying_type unit() { return underlying_type{ 1.,0. }; }
        static underlying_type append(underlying_type a, underlying_type b) {
            return underlying_type{ a.first * b.first, b.first * a.second + b.second };
        }
    };

    template <typename under = ll, under uni = 0, typename F = decltype(plus<ll>())>
    struct monoid_t {
        using underlying_type = under;
        static underlying_type unit() { return uni; }
        static underlying_type append(underlying_type a, underlying_type b) {
            return F(a, b);
        }
        static underlying_type act(underlying_type a, underlying_type b) {
            return F(a, b);
        }
    };


}

template<typename T>
struct AddAct :T {
    static typename T::underlying_type act(typename T::underlying_type a, typename T::underlying_type b) {
        return T::append(a, b);
    }
};

