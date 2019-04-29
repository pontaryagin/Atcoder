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

	template <typename T>
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

	template <typename T>
	struct min_t {
		typedef T underlying_type;
		static underlying_type unit() { return numeric_limits<T>::max(); }
		static underlying_type append(underlying_type a, underlying_type b) { return min(a, b); }
	};

	template <typename T>
	struct max_t {
		typedef T underlying_type;
		static underlying_type unit() { return numeric_limits<T>::min(); }
		static underlying_type append(underlying_type a, underlying_type b) { return max(a, b); }
	};

	struct linear_t {
		typedef pd underlying_type;
		static underlying_type unit() { return underlying_type{ 1.,0. }; }
		static underlying_type append(underlying_type a, underlying_type b) {
			return underlying_type{ a.first * b.first, b.first * a.second + b.second };
		}
	};


}


	


