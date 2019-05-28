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

	segment_tree(ll a_n) : size_original(a_n)
	{
		vector<underlying_type> initial_value = vector<underlying_type>(a_n, Monoid::unit());
		segment_tree_impl(a_n, initial_value);
	}

	segment_tree(ll a_n, vector<underlying_type>& initial_value) : size_original(a_n)
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
		underlying_type lacc = Monoid::unit(), racc = Monoid::unit();
		assert(l <= r && r <= n);
		l += n; r += n;
		for (; l < r; l /= 2, r /= 2) { // 1-based loop, 2x faster than recursion
			if (l % 2 == 1) lacc = Monoid::append(lacc, a[(l++) - 1]);
			if (r % 2 == 1) racc = Monoid::append(a[(--r) - 1], racc);
		}
		return Monoid::append(lacc, racc);
	}

	ll size() { return size_original; }

private:
	ll size_original;
	ll n;
	vector<underlying_type> a;
	void segment_tree_impl(ll a_n, vector<underlying_type>& initial_value)
	{
		assert(a_n == initial_value.size());
		n = 1; while (n < a_n) n *= 2;
		a.resize(2 * n - 1, Monoid::unit());
		rep(i, 0, initial_value.size()) {
			a[i + (n - 1)] = initial_value[i];
		}
		rrep(i, 0, n - 1) a[i] = Monoid::append(a[2 * i + 1], a[2 * i + 2]); // propagate initial values
	}


};


namespace Monoid{

	template <typename T = ll>
	struct min_indexed_t {
		typedef pair<T, ll> underlying_type;
		static underlying_type make_indexed(vector<T> v)
		{
			underlying_type w(v.size());
			rep(i, 0, v.size()) {
				w[i] = { v[i],i };
			}
			return w;
		}
		static underlying_type unit() { return make_pair(numeric_limits<T>::max(), -1); }
		static underlying_type append(underlying_type a, underlying_type b) { return min(a, b); }
	};

	template <typename T = ll>
	struct min_t {
		typedef T underlying_type;
		static underlying_type unit() { return numeric_limits<T>::max(); }
		static underlying_type append(underlying_type a, underlying_type b) { return min(a, b); }
	};

	template <typename T = ll>
	struct max_t {
		typedef T underlying_type;
		static underlying_type unit() { return numeric_limits<T>::min(); }
		static underlying_type append(underlying_type a, underlying_type b) { return max(a, b); }
	};

	template <typename T = ll>
	struct sum_t {
		typedef T underlying_type;
		static underlying_type unit() { return 0; }
		static underlying_type append(underlying_type a, underlying_type b) { return a+b; }
		static underlying_type iterate(underlying_type a, int n) { return a * n; }
	};

	struct linear_t {
		typedef pd underlying_type;
		static underlying_type unit() { return underlying_type{ 1.,0. }; }
		static underlying_type append(underlying_type a, underlying_type b) {
			return underlying_type{ a.first * b.first, b.first * a.second + b.second };
		}
	};


}

// 1) E is acting on T and 2) both should be monoid and 3) the action preserving monoid structure.
template <typename Monoid = void
	, typename T = typename Monoid::underlying_type
	, typename E = typename Monoid::underlying_type
	>
struct LazySegmentTree {
	int n;
	function<T(T, T)> f;
	function<T(T, E)> g;
	function<E(E, E)> h;
	function<E(E, int)> p;
	T t0;
	E e0;
	vector<T> dat;
	vector<E> laz;

	// Monoid has append, unit, iterate functions.
	//template<typename Monoid>  
	LazySegmentTree(int n_, vector<T> v = vector<T>())
		:f(Monoid::append), g(Monoid::append), h(Monoid::append), 
		t0(Monoid::unit()), e0(Monoid::unit()),p(Monoid::iterate)
	{
		init(n_);
		if (n_ == (int)v.size()) build(n_, v);
	}
	LazySegmentTree(int n_, function<T(T, T)> f, function<T(T, E)> g, function<E(E, E)> h, 
		T t0, E e0,	vector<T> v = vector<T>(), function<E(E, int)> p = [](E a, int) {return a;})
		:f(f), g(g), h(h), t0(t0), e0(e0), p(p) 
	{
		init(n_);
		if (n_ == (int)v.size()) build(n_, v);
	}
	void init(int n_) {
		n = 1;
		while (n < n_) n *= 2;
		dat.clear();
		dat.resize(2 * n - 1, t0);
		laz.clear();
		laz.resize(2 * n - 1, e0);
	}
	void build(int n_, vector<T> v) {
		for (int i = 0; i < n_; i++) dat[i + n - 1] = v[i];
		for (int i = n - 2; i >= 0; i--)
			dat[i] = f(dat[i * 2 + 1], dat[i * 2 + 2]);
	}
	inline void eval(int len, int k) {
		if (laz[k] == e0) return;
		if (k * 2 + 1 < n * 2 - 1) {
			laz[k * 2 + 1] = h(laz[k * 2 + 1], laz[k]);
			laz[k * 2 + 2] = h(laz[k * 2 + 2], laz[k]);
		}
		dat[k] = g(dat[k], p(laz[k], len));
		laz[k] = e0;
	}
	T update(int a, int b, E x, int k, int l, int r) {
		eval(r - l, k);
		if (r <= a || b <= l) return dat[k];
		if (a <= l && r <= b) {
			laz[k] = h(laz[k], x);
			return g(dat[k], p(laz[k], r - l));
		}
		return dat[k] = f(update(a, b, x, k * 2 + 1, l, (l + r) / 2),
			update(a, b, x, k * 2 + 2, (l + r) / 2, r));
	}
	T update(int a, int b, E x) {
		return update(a, b, x, 0, 0, n);
	}
	T query(int a, int b, int k, int l, int r) {
		eval(r - l, k);
		if (r <= a || b <= l) return t0;
		if (a <= l && r <= b) return dat[k];
		T vl = query(a, b, k * 2 + 1, l, (l + r) / 2);
		T vr = query(a, b, k * 2 + 2, (l + r) / 2, r);
		return f(vl, vr);
	}
	T query(int a, int b) {
		return query(a, b, 0, 0, n);
	}
};
	


