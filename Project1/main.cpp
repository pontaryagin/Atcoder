
//#pragma GCC optimize ("-O3")
#include <iostream>
#include <cmath>
#include <vector>
#include <stack>
#include <cstdio>
#include <string>
#include <bitset>
#include <list>
#include <set>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <functional>
#include <queue>
#include <regex>
#include <cassert>
#include <map>
#include <type_traits>
#include <array>
#include <cassert>
#include <typeinfo>
#include <time.h>
#include <iomanip>
#include <random>
#ifdef _MSC_VER
#include <intrin.h>
#define popcnt __popcnt64
//#  define __builtin_popcount __popcnt
#else
#define popcnt __builtin_popcountll
#endif
//#include "boost/variant.hpp"



using namespace std;

typedef long long ll;
constexpr ll MOD = 1000000007;
constexpr ll INF = 1LL << 60;

#define rep(i, N, M) for(ll i=N, i##_len=(M); i<i##_len; ++i)
#define rep_skip(i, N, M, ...) for(ll i=N, i##_len=(M); i<i##_len; i+=(skip))
#define rrep(i, N, M)  for(ll i=(M)-1, i##_len=(N-1); i>i##_len; --i)
#define pb push_back

typedef pair<double, double> pd;
typedef pair<ll, ll> pll;

template<int n>
struct tll_impl {
	using type = decltype(tuple_cat(tuple<ll>(), declval<typename tll_impl<n - 1>::type>()));
};
template<>
struct tll_impl<1> {
	using type = tuple<ll>;
};
template<int n>
using tll = typename tll_impl<n>::type;

template<class T>
constexpr ll SZ(T& v) { return static_cast<ll>(v.size()); };

template<int n, typename T>
struct vec_t_impl {
	using type = vector<typename vec_t_impl<n-1,T>::type>;
};
template<typename T>
struct vec_t_impl<1,T> {
	using type = vector<T>;
};
template<int n, typename T>
using vec_t = typename vec_t_impl<n, T>::type;
// check 
static_assert(is_same<vec_t<3,ll>, vector<vector<vector<ll>>>>::value, "");

// decompose vector into basetype and dimension.
template<typename T> 
struct vec_dec {
	static constexpr int dim = 0;
	using type  = T;
};
template<typename T>
struct vec_dec<vector<T>> {
	static constexpr int dim = vec_dec<T>::dim+1;
	using type  = typename vec_dec<T>::type;
};
static_assert(is_same<typename vec_dec<vec_t<3, ll>>::type, ll>::value, "");
static_assert(vec_dec<vec_t<3, ll>>::dim == 3, "");

template<typename T = ll>
vector<T> make_v(size_t a) { return vector<T>(a); }

template<typename T = ll, typename... Ts>
auto make_v(size_t a, Ts... ts) {
	return vector<decltype(make_v<T>(ts...))>(a, make_v<T>(ts...));
}
// ex:  auto dp =  make_v<ll>(4,5) => vector<vector<ll>> dp(4,vector<ll>(5));

// check if T is vector
template < typename T >
struct is_vector : std::false_type {};

template < typename T >
struct is_vector<vector<T>> : std::true_type {};
static_assert(is_vector<vector<ll>>::value == true && is_vector<ll>::value == false, "");

// check if T is vector
template < typename T>
struct is_pair : std::false_type {};

template < typename T, typename S >
struct is_pair<pair<T, S>> : std::true_type {};
static_assert(is_pair<pll>::value == true && is_pair<ll>::value == false, "");

template<typename T, typename V, typename enable_if<!is_vector<T>::value, nullptr_t>::type = nullptr>
void fill_v(T& t, const V& v) { t = v; }

template<typename T, typename V, typename enable_if<is_vector<T>::value, nullptr_t>::type = nullptr>
void fill_v(T& t, const V& v) {
	for (auto &&x : t)
		fill_v(x, v);
}
// ex:  fill_v(dp, INF);

template<typename T, typename enable_if < !is_vector<T>::value && !is_pair<T>::value, nullptr_t > ::type = nullptr >
void read_v(T& x) {	cin >> x;}

template<typename T, typename enable_if<is_pair<T>::value, nullptr_t>::type = nullptr>
void read_v(T& x) { read_v(x.first); read_v(x.second); }

template<typename T, typename enable_if<is_vector<T>::value, nullptr_t>::type = nullptr>
void read_v(T& x) { rep(i,0,x.size()) read_v(x[i]); }

template<typename T, typename Delim_t = string, typename enable_if<!is_vector<T>::value, nullptr_t>::type = nullptr>
void write_v(T & x, Delim_t delim = " ") { cout << x << delim; }

template<typename T, typename Delim_t = string, typename enable_if<is_vector<T>::value, nullptr_t>::type = nullptr>
void write_v(T& x, Delim_t delim = " ") { rep(i, 0, x.size()) write_v(x[i]); cout << endl; }

typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<pll> vpll;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef vector<string> vs;
template<typename T>
using pq_greater = priority_queue<T, vector<T>, greater<T>>;


#define all(a)  (a).begin(),(a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define perm(c) sort(all(c));for(bool c##perm=1;c##perm;c##perm=next_permutation(all(c)))

template<typename T> void chmin(T &a, T b) {
	if (a > b) a = b;
}
template<typename T> void chmax(T &a, T b) {
	if (a < b) a = b;
}

vll seq(ll i, ll j) {
	vll res(j - i);
	rep(k, i, j) res[k] = i + k;
	return res;
}

constexpr ll POW_0(ll x, ll y) {
	if (y == 0)return 1;
	if (y == 1)return x ;
	if (y == 2)return x * x ;
	if (y % 2 == 0)return POW_0(POW_0(x, y / 2), 2LL);
	return ((POW_0(POW_0(x, y / 2), 2LL)) * (x)) ;
}

constexpr ll POW(ll x, ll y, ll mod = MOD) {
	if (mod == 0)return POW_0(x, y);
	if (y == 0)return 1;
	if (y == 1)return x % mod;
	if (y == 2)return x * x % mod;
	if (y % 2 == 0)return POW(POW(x, y / 2, mod), 2LL, mod) % mod;
	return ((POW(POW(x, y / 2, mod), 2LL, mod)) * (x % mod)) % mod;
}

template<
	typename Inputs,
	typename Functor,
	typename T = typename Inputs::value_type>
	void sort_by(Inputs& inputs, Functor f) {
	std::sort(std::begin(inputs), std::end(inputs),
		[&f](const T& lhs, const T& rhs) { return f(lhs) < f(rhs); });
}

template<
	typename Inputs,
	typename Functor,
	typename T = typename Inputs::value_type>
	void stable_sort_by(Inputs& inputs, Functor f) {
	std::stable_sort(std::begin(inputs), std::end(inputs),
		[&f](const T& lhs, const T& rhs) { return f(lhs) < f(rhs); });
}






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
	, typename E = typename Monoid::underlying_type>
struct LazySegmentTree {
	typedef function<T(T, T)> F;
	typedef function<T(T, E)> G;
	typedef function<E(E, E)> H;
	typedef function<E(E, int)> P;
	int n;
	F f;
	G g;
	H h;
	P p;
	T d1;
	E d0;
	vector<T> dat;
	vector<E> laz;

	// Monoid has append, unit, iterate functions.
	//template<typename Monoid>  
	LazySegmentTree(int n_, vector<T> v = vector<T>())
		:f(Monoid::append), g(Monoid::append), h(Monoid::append), 
		d1(Monoid::unit()), d0(Monoid::unit()),p(Monoid::iterate)
	{
		init(n_);
		if (n_ == (int)v.size()) build(n_, v);
	}
	LazySegmentTree(int n_, F f, G g, H h, T d1, E d0,
		vector<T> v = vector<T>(), P p = [](E a, int b) {return a; }) 
		:f(f), g(g), h(h), d1(d1), d0(d0), p(p) 
	{
		init(n_);
		if (n_ == (int)v.size()) build(n_, v);
	}
	void init(int n_) {
		n = 1;
		while (n < n_) n *= 2;
		dat.clear();
		dat.resize(2 * n - 1, d1);
		laz.clear();
		laz.resize(2 * n - 1, d0);
	}
	void build(int n_, vector<T> v) {
		for (int i = 0; i < n_; i++) dat[i + n - 1] = v[i];
		for (int i = n - 2; i >= 0; i--)
			dat[i] = f(dat[i * 2 + 1], dat[i * 2 + 2]);
	}
	inline void eval(int len, int k) {
		if (laz[k] == d0) return;
		if (k * 2 + 1 < n * 2 - 1) {
			laz[k * 2 + 1] = h(laz[k * 2 + 1], laz[k]);
			laz[k * 2 + 2] = h(laz[k * 2 + 2], laz[k]);
		}
		dat[k] = g(dat[k], p(laz[k], len));
		laz[k] = d0;
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
		if (r <= a || b <= l) return d1;
		if (a <= l && r <= b) return dat[k];
		T vl = query(a, b, k * 2 + 1, l, (l + r) / 2);
		T vr = query(a, b, k * 2 + 2, (l + r) / 2, r);
		return f(vl, vr);
	}
	T query(int a, int b) {
		return query(a, b, 0, 0, n);
	}
};
	















// ============================ Header  =================================

template<typename Inputs>
void sort_uniq(Inputs& inputs) {
	sort(all(inputs));
	inputs.erase(unique(all(inputs)), inputs.end());
}

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	
	int n, q;
	cin >> n >> q;
	LazySegmentTree<Monoid::sum_t<>> ch(n);
	for (int i = 0; i < q; i++) {
		int c, s, t, x;
		cin >> c;
		if (c) {
			cin >> s >> t;
			cout << ch.query(s - 1, t) << endl;
		}
		else {
			cin >> s >> t >> x;
			ch.update(s - 1, t, x);
		}
	}

	
	return 0;

}

