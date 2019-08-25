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

	underlying_type query(ll i) { // return value at i
		assert(0 <= i && i < size_original);
		return a[i + n - 1];
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


namespace M {

	template <typename T = ll>
	struct sum_t {
		typedef T underlying_type;
		static underlying_type unit() { return 0; }
		static underlying_type append(underlying_type a, underlying_type b) { return a + b; }
		static underlying_type iterate(underlying_type a, int n) { return a * n; }
	};

	template<typename S, typename T>
	struct pair_t {
		typedef pair<typename S::underlying_type, typename T::underlying_type> underlying_type;
		static underlying_type unit() { return make_pair(S::unit(), T::unit()); }
		static underlying_type append(underlying_type a, underlying_type b) { return make_pair(S::append(a.first, b.first), T::append(a.second, b.second)); }
		static underlying_type iterate(underlying_type a, int n) { return make_pair(S::iterate(a.first,n), T::iterate(a.second,n)); }
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

// 1) E is acting on T and 2) both should be monoid and 3) the action preserving monoid structure.
// requires 
template <typename Monoid, typename ActionMonoid = AddAct<Monoid>>
struct LazySegmentTree {
	int n;
	using M = typename Monoid::underlying_type;
	using E = typename ActionMonoid::underlying_type;
	function<M(M, M)> f = Monoid::append;
	function<M(M, E)> act = ActionMonoid::act;

	function<E(E, E)> h = ActionMonoid::append;
	function<E(E, int)> iterate = ActionMonoid::iterate;
	M m0 = Monoid::unit();
	E e0 = ActionMonoid::unit();
	vector<M> dat;
	vector<E> laz;

	// Monoid has append, unit, iterate functions.
	//template<typename Monoid>  
	LazySegmentTree(int n_, vector<M> v = vector<M>())
	{
		init(n_);
		if (n_ == (int)v.size()) build(n_, v);
	}

	LazySegmentTree(int n_, function<M(M, M)> f, function<M(M, E)> act, function<E(E, E)> h,
		M m0, E e0, vector<M> v = vector<M>(), function<E(E, int)> iterate = [](E a, int) {return a; })
		:f(f), act(act), h(h), m0(m0), e0(e0), iterate(iterate)
	{
		init(n_);
		if (n_ == (int)v.size()) build(n_, v);
	}
	void init(int n_) {
		n = 1;
		while (n < n_) n *= 2;
		dat.clear();
		dat.resize(2 * n - 1, m0);
		laz.clear();
		laz.resize(2 * n - 1, e0);
	}
	void build(int n_, vector<M> v) {
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
		dat[k] = act(dat[k], iterate(laz[k], len));
		laz[k] = e0;
	}
	M update(int a, int b, E x, int k, int l, int r) {
		eval(r - l, k);
		if (r <= a || b <= l) return dat[k];
		if (a <= l && r <= b) {
			laz[k] = h(laz[k], x);
			return act(dat[k], iterate(laz[k], r - l));
		}
		return dat[k] = f(update(a, b, x, k * 2 + 1, l, (l + r) / 2),
			update(a, b, x, k * 2 + 2, (l + r) / 2, r));
	}
	M update(int a, int b, E x) {
		return update(a, b, x, 0, 0, n);
	}
	M update(int a,  E x) {
		return update(a, a + 1, x);
	}
	M query(int a, int b, int k, int l, int r) {
		eval(r - l, k);
		if (r <= a || b <= l) return m0;
		if (a <= l && r <= b) return dat[k];
		M vl = query(a, b, k * 2 + 1, l, (l + r) / 2);
		M vr = query(a, b, k * 2 + 2, (l + r) / 2, r);
		return f(vl, vr);
	}
	M query(int a, int b) {
		return query(a, b, 0, 0, n);
	}
};



