
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
#include <sstream>
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
#define fir first
#define sec second

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
vector<T> makev(size_t a) { return vector<T>(a); }

template<typename T = ll, typename... Ts>
auto makev(size_t a, Ts... ts) {
	return vector<decltype(makev<T>(ts...))>(a, makev<T>(ts...));
}
// ex:  auto dp =  makev<ll>(4,5) => vector<vector<ll>> dp(4,vector<ll>(5));

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
void read(T& x) {	cin >> x;}

template<typename T, typename enable_if<is_pair<T>::value, nullptr_t>::type = nullptr>
void read(T& x) { read(x.first); read(x.second); }

template<typename T, typename enable_if<is_vector<T>::value, nullptr_t>::type = nullptr>
void read(T& x) { rep(i,0,x.size()) read(x[i]); }

template<typename T, typename Delim_t = string, typename enable_if<!is_vector<T>::value, nullptr_t>::type = nullptr>
void write(T & x, Delim_t delim = " ") { cout << x << delim; }

template<typename T, typename Delim_t = string, typename enable_if<is_vector<T>::value, nullptr_t>::type = nullptr>
void write(T& x, Delim_t delim = " ") { rep(i, 0, x.size()) write(x[i], (i == (x.size() - 1) ? "" : delim)); cout << '\n'; }

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

constexpr ll POW(ll x, ll y, ll mod = 0) {
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

template<typename Inputs>
void sort_uniq(Inputs& inputs) {
	sort(all(inputs));
	inputs.erase(unique(all(inputs)), inputs.end());
}

vector<string> split(const string& s, char delim) {
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		if (!item.empty()) {
			elems.push_back(item);
		}
	}
	return elems;
}







ll div_ferm(ll val, ll  b, ll mod) {
	return (val* POW(b, mod - 2, mod)) % mod;
}


// === Modint ===
static uint_fast64_t runtime_modulus = MOD;

template <std::uint_fast64_t modulus = MOD> 
class modint 
{
	using u64 = std::uint_fast64_t;

public:
	u64 val; 

	template<ll modulus_ = modulus, enable_if_t<(modulus_ <= 0), nullptr_t> = nullptr >
	auto mod() {  return runtime_modulus; }
	template<ll modulus_ = modulus, enable_if_t<(modulus_ > 0), nullptr_t> = nullptr >
	auto mod() { return modulus; }
	//template<class Ret = u64 &>
	//static auto modulo() -> std::enable_if_t<(modulus <= 0), Ret> { static u64 runtime_modulus=0; return runtime_modulus; }
	//template<class Ret = const u64>
	//static auto mod() -> std::enable_if_t<(modulus <= 0), Ret> { return modulo(); }
	//template<class Ret = const u64>
	//static constexpr auto mod()->std::enable_if_t<(modulus > 0), Ret> { return modulus; }
	//template<ll modulus_ = modulus, enable_if_t<(modulus_ <= 0), nullptr_t> = nullptr >
	//static void set_modulo(u64 mod) { runtime_modulus = mod; }

	constexpr modint(): val(0){}
	constexpr modint(const u64 x) noexcept : val(x % mod()) {}
	constexpr modint inv() { return pow(mod() - 2); }

	constexpr u64 value() const noexcept { return val; }
	constexpr modint operator+(const modint rhs) const noexcept {
		return modint(*this) += rhs;
	}
	constexpr modint operator-(const modint rhs) const noexcept {
		return modint(*this) -= rhs;
	}
	constexpr modint operator*(const modint rhs) const noexcept {
		return modint(*this) *= rhs;
	}
	constexpr modint operator/(const modint rhs) const noexcept {
		return modint(*this) /= rhs;
	}
	modint &operator+=(const modint rhs) noexcept {
		val += rhs.val;
		if (val >= mod()) {
			val -= mod();
		}
		return *this;
	}
	modint &operator-=(const modint rhs) noexcept {
		if (val < rhs.val) {
			val += mod();
		}
		val -= rhs.val;
		return *this;
	}
	modint &operator*=(const modint rhs) noexcept {
		val = val * rhs.val % mod();
		return *this;
	}
	modint &operator/=(modint rhs) noexcept {
		u64 exp = mod() - 2;
		while (exp) {
			if (exp % 2) {
				*this *= rhs;
			}
			rhs *= rhs;
			exp /= 2;
		}
		return *this;
	}
	modint &operator++() noexcept {
		return *this += modint(1);
	}
	modint &operator++(int) noexcept {
		auto t = *this;
		*this += modint(1);
		return t;
	}
	modint &operator--() noexcept {
		return *this -= modint(1);
	}
	modint &operator--(int) noexcept {
		auto t = *this;
		*this -= modint(1);
		return t;
	}
	constexpr modint operator-() { return val ? mod() - val : val; }
	constexpr bool operator==(const modint rhs) const noexcept { return val == rhs.value(); }
	constexpr bool operator!=(const modint rhs)const  noexcept { return val != rhs.value(); }
	constexpr bool operator <(const modint rhs)const  noexcept { return val < rhs.value(); }
	static constexpr modint zero() { return modint(0); }
	static constexpr modint unit() { return modint(1); }

	modint pow(long long k) const {
		modint v = *this;
		modint res(1), tmp(v);
		while (k) {
			if (k & 1) res *= tmp;
			tmp *= tmp;
			k >>= 1;
		}
		return res;
	}
	
	u64 log(modint b) {
		modint val = *this;
		const u64 sq = 40000;
		map<modint, u64> dp;
		//dp.reserve(sq);
		modint res(1);
		for (ll r = 0; r < sq; r++) {
			if (!dp.count(res)) dp[res] = r;
			res *= val;
		}
		modint p = val.inv().pow(sq);
		res = b;
		for (ll q = 0; q <= mod() / sq + 1; q++) {
			if (dp.count(res)) {
				u64 idx = q * sq + dp[res];
				if (idx > 0) return idx;
			}
			res *= p;
		}
		return INF;
	}
	friend ostream& operator <<(ostream& o, const modint<modulus>& t) {
		o << t.value();
		return o;
	}
	friend istream& operator >>(istream& in, modint<modulus>& t) {
		uint_fast64_t x;
		in >> x;
		t = modint<modulus>(x);
		return in;
	}
	friend modint<modulus> POW(modint<modulus> x, ll n) {
		return modint<modulus>(POW(x.value(), n, mod()));
	}

};


class Combination {
	// this calculates combination (nCk).
	// Constructor runs in O(MAX).
	// get(n,k) returns nCk in O(1).

	ll N_MAX, mod;
	vll fac;
	vll finv;
	vll inv;
public:
	Combination(ll mod = MOD, ll N_MAX = 210000)
		:mod(mod), N_MAX(max(N_MAX, 2LL)), fac(vll(N_MAX + 1)), finv(vll(N_MAX + 1)), inv(vll(N_MAX + 1)) {
		fac[0] = fac[1] = 1;
		finv[0] = finv[1] = 1;
		inv[1] = 1;
		pre_process(2LL, N_MAX + 1);
	}

	ll operator()(ll n, ll k) {
		// choose k from n
		if (N_MAX < n)
			pre_process(N_MAX + 1, n + 1);

		if (n < k)return 0;
		if (n < 0 || k < 0)return 0;
		return fac[n] * (finv[k] * finv[n - k] % mod) % mod;
	}
	ll H(ll n, ll k) {
		// 1) 区間[0, k) を（空を許して）n個に分割する場合の数
		// 2) n個の中からk個を重複を許して選ぶ
		return (n==0 && k==0)? 1 : operator()(n + k - 1, k);
	}
private:
	void pre_process(ll m, ll n) {
		if (N_MAX < n) {
			fac.resize(n); inv.resize(n); finv.resize(n);
		}
		rep(i, m, n) {
			fac[i] = fac[i - 1] * i % mod;
			inv[i] = mod - inv[mod%i] * (mod / i) % mod;
			finv[i] = finv[i - 1] * inv[i] % mod;
		}
	}
};


ll choose(int n, int r) { // O(r) for small n
	ll acc = 1;
	rep(i, 0, r) acc = acc * (n - i) / (i + 1);
	return acc;
}

ll gcd(ll val, ll b) {
	if (val%b == 0) return b;
	else return gcd(b, val%b);
}




vll getDivisors(ll n) {
	vll res;
	ll i = 1;

	for (; i*i < n; i++) {
		if (n%i == 0) {
			res.push_back(i);
			res.push_back(n / i);
		}
	}
	if (i*i == n)res.push_back(i);
	sort(res.begin(), res.end());
	return res;
}

vll getDivisors(ll n, ll m) {
	// O(sqrt(min(n,m)))
	if (n > m) swap(n, m);
	vll res;
	ll i = 1;

	for (; i*i < n; i++) {
		if (n%i == 0) {
			if (m%i == 0) res.push_back(i);
			if (m % (n / i) == 0) res.push_back(n / i);
		}
	}
	if (i*i == n) if (m%i == 0) res.push_back(i);
	sort(res.begin(), res.end());
	return res;
}

vector<pll >prime_factorize(ll n) {
	vector<pll> res;
	for (ll p = 2; p*p <= n; ++p) {
		if (n%p != 0)continue;
		ll num = 0;
		while (n%p == 0) { ++num; n /= p; }
		res.push_back({ p,num });
	}
	if (n != 1) res.push_back(make_pair(n, 1));
	return res;
}














// ============================ Header  =================================


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);
	
	ll n; cin >> n;
	vll a(n);
	read(a);
	vll res(n);
	vll S0(n + 1), S1(n+1);
	rep(i, 0, n / 2) {
		S0[i*2+2 ] =S0[i*2]+ a[i * 2];
	}
	rep(i, 0, n / 2) {
		S1[1+ i * 2 + 2] = S1[1+ i * 2] + a[1+ i * 2];
	}

	ll tot = accumulate(all(a), 0LL);
	res[n - 1] = tot;
	rep(i, 0, n / 2)
		res[n - 1] -= a[i * 2]*2;
	ll x = n - 1;
	rep(i, 0, n-1) {
		ll tmp = tot ;
		ll y = (x - 1 >= 0 ? x - 1 : x + n - 1);

		tmp -= a[y]*2;
		tmp = tot - tmp;
		tmp -= res[x];
		res[y] = tmp;
		x = y;



		/*if (i == 0) {
			res[0] = tot -S1[n];
		}
		else if(i==1){
			res[1] = tot;
			res[1] -= a[0];
			res[1] -= S1[n];
		}
		else if(i%2 ==0){
			res[i] = tot;
			res[i] -= S1[i];
			res[i] -= (S1
		}*/
	}
	cout << res[0];
	rep(i, 1, n ) {
		cout << " " << res[i] ;
	}
	cout << endl;
	return 0;

}
