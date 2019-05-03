

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




using namespace std;

typedef long long ll;
constexpr ll MOD = 1000000007;
constexpr ll INF = 1LL << 62;

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

template<typename T, typename V>
typename enable_if<!is_vector<T>::value>::type
fill_v(T& t, const V& v) { t = v; }

template<typename T, typename V>
typename enable_if<is_vector<T>::value>::type
fill_v(T& t, const V& v) {
	for (auto &&x : t)
		fill_v(x, v);
}
// ex:  fill_v(dp, INF);

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

constexpr ll POW_(ll n, ll m) {
	ll res = 1;
	rep(i, 0, m) {
		res *= n;
	}
	return res;
}

template<ll mod = 0>
constexpr ll POW(ll x, ll n) {
	if (x == 2)
	{
		return (1LL << n) % mod;
	}
	if (n == 0)return 1;
	if (n == 1)return x % mod;
	if (n % 2 == 0)return POW_(POW<mod>(x, n / 2), 2LL) % mod;
	return ((POW_(POW<mod>(x, n / 2), 2LL) % mod)*(x%mod)) % mod;
}
template<>
constexpr ll POW<0>(ll x, ll n) {
	if (x == 2)
	{
		return 1LL << n;
	}
	if (n == 0)return 1;
	if (n == 1)return x;
	if (n % 2 == 0) return POW_(POW(x, n / 2), 2);
	return (POW_(POW(x, n / 2), 2))*x;
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





ll POW(ll n, ll m) {
	ll res = 1;
	rep(i, 0, m) {
		res *= n;
	}
	return res;
}

ll bin_power(ll x, ll y, ll mod) {
	if (y == 0)return 1;
	if (y == 1)return x % mod;
	if (y % 2 == 0)return POW(bin_power(x, y / 2, mod), 2LL) % mod;
	return ((POW(bin_power(x, y / 2, mod), 2LL) % mod)*(x%mod)) % mod;
}
ll div_ferm(ll a, ll  b, ll mod) {
	return (a* bin_power(b, mod - 2, mod)) % mod;
}


// === Modint ===

template <std::uint_fast64_t Modulus = 1000000007> 
class modint 
{
	using u64 = std::uint_fast64_t;

public:
	u64 a;

	constexpr modint(const u64 x = 0) noexcept : a(x % Modulus) {}
	//constexpr modint(const modint& rhs) noexcept {
	//	this->a = rhs.value();
	//}
	//constexpr modint &operator=(const modint &rhs) noexcept {
	//	this->a = rhs.value();
	//	return *this;
	//}
	constexpr u64 value() const noexcept { return a; }
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
		a += rhs.a;
		if (a >= Modulus) {
			a -= Modulus;
		}
		return *this;
	}
	modint &operator-=(const modint rhs) noexcept {
		if (a < rhs.a) {
			a += Modulus;
		}
		a -= rhs.a;
		return *this;
	}
	modint &operator*=(const modint rhs) noexcept {
		a = a * rhs.a % Modulus;
		return *this;
	}
	modint &operator/=(modint rhs) noexcept {
		u64 exp = Modulus - 2;
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

};
template<uint_fast64_t Modulus>
ostream& operator <<(ostream &o, const modint<Modulus> &t) {
	o << t.value();
	return o;
}
template<uint_fast64_t Modulus>
istream& operator >>(istream &in, modint<Modulus> &t) {
	uint_fast64_t x;
	in >> x;
	t = modint<Modulus>(x);
	return in;
}
template<uint_fast64_t Modulus>
modint<Modulus> POW(modint<Modulus> x, ll n) {
	return modint<Modulus>(POW<Modulus>(x.value(), n));
}

// === Mll ===

template<typename T, T MOD = 1000000007>
struct Mint {
	T v;
	Mint() :v(0) {}
	Mint(signed v) :v(v) {}
	Mint(long long t) { v = t % MOD; if (v < 0) v += MOD; }




	Mint inv() { return pow(MOD - 2); }

	Mint& operator+=(Mint a) { v += a.v; if (v >= MOD)v -= MOD; return *this; }
	Mint& operator-=(Mint a) { v += MOD - a.v; if (v >= MOD)v -= MOD; return *this; }
	Mint& operator*=(Mint a) { v = 1LL * v*a.v%MOD; return *this; }
	Mint& operator/=(Mint a) { return (*this) *= a.inv(); }

	Mint operator+(Mint a) const { return Mint(v) += a; };
	Mint operator-(Mint a) const { return Mint(v) -= a; };
	Mint operator*(Mint a) const { return Mint(v) *= a; };
	Mint operator/(Mint a) const { return Mint(v) /= a; };

	Mint operator-() { return v ? MOD - v : v; }

	bool operator==(const Mint a)const { return v == a.v; }
	bool operator!=(const Mint a)const { return v != a.v; }
	bool operator <(const Mint a)const { return v < a.v; }

	static Mint pow(Mint v, long long k) {
		Mint res(1), tmp(v);
		while (k) {
			if (k & 1) res *= tmp;
			tmp *= tmp;
			k >>= 1;
		}
		return res;
	}

	// find x s.t. a^x = b
	static T log(Mint a, Mint b) {
		const T sq = 40000;
		unordered_map<T, T> dp;
		dp.reserve(sq);
		Mint res(1);
		for (ll r = 0; r < sq; r++) {
			if (!dp.count(res)) dp[res] = r;
			res *= a;
		}
		Mint p = pow(a.inv(), sq);
		res = b;
		for (ll q = 0; q <= MOD / sq + 1; q++) {
			if (dp.count(res)) {
				T idx = q * sq + dp[res];
				if (idx > 0) return idx;
			}
			res *= p;
		}
		return T(-1);
	}

	static vector<Mint> fact, finv, invs;

	static void init(ll n) {
		if (n + 1 <= (signed)fact.size()) return;
		fact.assign(n + 1, 1);
		finv.assign(n + 1, 1);
		invs.assign(n + 1, 1);

		for (ll i = 1; i <= n; i++) fact[i] = fact[i - 1] * Mint(i);
		finv[n] = Mint(1) / fact[n];
		for (ll i = n; i >= 1; i--) finv[i - 1] = finv[i] * Mint(i);
		for (ll i = 1; i <= n; i++) invs[i] = finv[i] * fact[i - 1];
	}

	static Mint comb(long long n, ll k) {
		Mint res(1);
		for (ll i = 0; i < k; i++) {
			res *= Mint(n - i);
			res /= Mint(i + 1);
		}
		return res;
	}

	static Mint C(ll n, ll k) {
		if (n < k || k < 0) return Mint(0);
		init(n);
		return fact[n] * finv[n - k] * finv[k];
	}

	static Mint P(ll n, ll k) {
		if (n < k || k < 0) return Mint(0);
		init(n);
		return fact[n] * finv[n - k];
	}

	static Mint H(ll n, ll k) {
		if (n < 0 || k < 0) return Mint(0);
		if (!n && !k) return Mint(1);
		init(n + k - 1);
		return C(n + k - 1, k);
	}

	static Mint S(ll n, ll k) {
		Mint res;
		init(k);
		for (ll i = 1; i <= k; i++) {
			Mint tmp = C(k, i)*Mint(i).pow(n);
			if ((k - i) & 1) res -= tmp;
			else res += tmp;
		}
		return res *= finv[k];
	}

	static vector<vector<Mint> > D(ll n, ll m) {
		vector<vector<Mint> > dp(n + 1, vector<Mint>(m + 1, 0));
		dp[0][0] = Mint(1);
		for (ll i = 0; i <= n; i++) {
			for (ll j = 1; j <= m; j++) {
				if (i - j >= 0) dp[i][j] = dp[i][j - 1] + dp[i - j][j];
				else dp[i][j] = dp[i][j - 1];
			}
		}
		return dp;
	}

	static Mint B(ll n, ll k) {
		Mint res;
		for (ll j = 1; j <= k; j++) res += S(n, j);
		return res;
	}

	static Mint montmort(ll n) {
		Mint res;
		init(n);
		for (ll k = 2; k <= n; k++) {
			if (k & 1) res -= finv[k];
			else res += finv[k];
		}
		return res *= fact[n];
	}

	static Mint LagrangePolynomial(vector<Mint> &y, Mint t) {
		ll n = y.size() - 1;
		if (t.v <= n) return y[t.v];
		init(n + 1);
		Mint num(1);
		for (ll i = 0; i <= n; i++) num *= t - Mint(i);
		Mint res;
		for (ll i = 0; i <= n; i++) {
			Mint tmp = y[i] * num / (t - Mint(i))*finv[i] * finv[n - i];
			if ((n - i) & 1) res -= tmp;
			else res += tmp;
		}
		return res;
	}
};


class Combination {
	// this calculates combination (nCk).
	// Constructor runs in O(MAX).
	// get(n,k) returns nCk in O(1).

	ll MAX, MOD;
	vll fac;
	vll finv;
	vll inv;
public:
	Combination(ll MAX = 210000, ll MOD = 1000000007)
		:MOD(MOD), MAX(max(MAX, 2LL)), fac(vll(MAX + 1)), finv(vll(MAX + 1)), inv(vll(MAX + 1)) {
		fac[0] = fac[1] = 1;
		finv[0] = finv[1] = 1;
		inv[1] = 1;
		pre_process(2LL, MAX + 1);
	}

	ll get(ll n, ll k) {
		if (MAX < n)
			pre_process(MAX + 1, n + 1);

		if (n < k)return 0;
		if (n < 0 || k < 0)return 0;
		return fac[n] * (finv[k] * finv[n - k] % MOD) % MOD;
	}
private:
	void pre_process(ll m, ll n) {
		if (MAX < n) {
			fac.resize(n); inv.resize(n); finv.resize(n);
		}
		rep(i, m, n) {
			fac[i] = fac[i - 1] * i % MOD;
			inv[i] = MOD - inv[MOD%i] * (MOD / i) % MOD;
			finv[i] = finv[i - 1] * inv[i] % MOD;
		}
	}
};


ll choose(int n, int r) { // O(r) for small n
	ll acc = 1;
	rep(i, 0, r) acc = acc * (n - i) / (i + 1);
	return acc;
}

ll gcd(ll a, ll b) {
	if (a%b == 0) return b;
	else return gcd(b, a%b);
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








struct Edge
{
	ll from;
	ll to;
	ll cost=1;
	Edge reverse() const {
		return Edge{ to, from , cost };
	}
	Edge(ll from , ll to, ll cost=1) : from(from),to(to),cost(cost){};
	Edge(pll e) { from = e.first; to = e.second; cost = 1; }
	Edge() :from(0), to(0), cost(0){ };
};

struct Graph {
	vector<Edge> edges;
	vector<vector<ll>> out_edges;
	vector<vector<ll>> in_edges;
	Graph(ll nodeSize, const vector<Edge>& edges = vector<Edge>()): out_edges(nodeSize), in_edges(nodeSize), edges(move(edges)){
		rep(i, 0, edges.size()) {
			in_edges[edges[i].to].push_back(i);
			out_edges[edges[i].from].push_back(i);
		}
	}
	Graph(ll nodeSize, vector<pll> edges_) : out_edges(nodeSize), in_edges(nodeSize), edges(edges_.size()) {
		rep(i, 0, edges.size()) {
			edges[i] = edges_[i];
			in_edges[edges[i].to].push_back(i);
			out_edges[edges[i].from].push_back(i);
		}
	}
	
	Edge& operator[](ll ind) {
		return this->edges[ind];
	}
	const Edge& operator[](ll ind) const {
		return this->edges[ind];
	}
	size_t size() const { return out_edges.size(); }
	void push(Edge edge){
		assert(max(edge.from, edge.to) < out_edges.size());
		edges.push_back(edge);
		out_edges[edge.from].push_back(edges.size() - 1);
		in_edges[edge.to].push_back(edges.size() - 1);
	}
	void erase(ll from, ll to) {
		// O(nodeSize)
		assert(max(from, to) < out_edges.size());
		for (ll k : out_edges[from]) {
			if (edges[k].to = to) {
				edges.erase(edges.begin() + k);
			}
		}
	}
	void erase(ll edge_ind) {
		// O(nodeSize)
		assert(edge_ind < edges.size());
		edges.erase(edges.begin() + edge_ind);
		auto& esin = in_edges[edges[edge_ind].to];
		esin.erase(find(all(esin), edge_ind));
		auto& esout = out_edges[edges[edge_ind].from];
		esout.erase(find(all(esout), edge_ind));
	}
	void push(vector<Edge> edges) {
		for (const Edge& e : edges) {
			push(e);
		}
	}

	vll get_topologically_sorted_nodes()
	{
		// graph needs to be represented by adjacent list.
		// complexity: O( node size + edge size)
		ll nodeSize = this->size();

		// find root
		vll roots;
		vll inDegree(nodeSize);
		rep(i, 0, nodeSize)
		{
			for (ll sibling_ind : this->out_edges[i]) {
				inDegree[this->edges[sibling_ind].to]++;
			}
		}


		rep(i, 0, nodeSize) {
			if (inDegree[i] == 0) {
				roots.push_back(i);
			}
		}

		stack<ll> parents;
		for (ll i : roots)
			parents.push(i);

		vll sortedNodes;
		while (!parents.empty()) {
			ll parent = parents.top();
			parents.pop();
			sortedNodes.push_back(parent);
			for (ll sibling_ind : this->out_edges[parent]) {
				Edge sibling = this->edges[sibling_ind];
				inDegree[sibling.to]--;
				if (inDegree[sibling.to] == 0) {
					parents.push(sibling.to);
				}
			}
		}
		return sortedNodes;
	}
	void topological_sort() {
		vll sorted = get_topologically_sorted_nodes();
		vll new_ind(sorted.size());
		vector<Edge> new_edges;
		rep(i, 0, sorted.size()) {
			new_ind[sorted[i]] = i;
		}
		for (Edge& e : edges) {
			new_edges.emplace_back(Edge{ new_ind[e.from], new_ind[e.to],e.cost });
		}
		*this = Graph(this->size(), new_edges);
	}

};

pair<vll, vll> dijkstra(const Graph& graph, size_t start) {
	// graph: weighted directed graph of adjacent representation
	// start: index of start point
	// return1: minimum path length from start
	// return2: concrete shortest path info
	ll node_size = graph.size();
	vll dist(node_size, 1LL << 60);
	vll from_list(node_size, -1);
	dist[start] = 0;
	pq_greater<pair<ll, pll>> pq;
	pq.push({ 0, {start, start} });
	while (!pq.empty()) {
		auto node = pq.top(); pq.pop();
		// if not shortest path fixed, fix
		ll from = node.second.first;
		ll to = node.second.second;
		if (from_list[to] != -1)
			continue;
		from_list[to] = from;

		for (ll edge_ind : graph.out_edges[to]) {
			const Edge& edge = graph[edge_ind];
			ll adj = edge.to;
			ll cost = dist[to] + edge.cost;
			if (dist[adj] > cost) {
				dist[adj] = min(dist[adj], cost);
				pq.push({ cost ,{to, adj} });
			}
		}
	}
	return { dist, from_list };

}

vll shortest_path(const vll& from_list, ll start, ll goal) {
	// usage : vll path =  shortest_path(dijkstra(g,s).second, s, g);
	vll path;
	path.emplace_back(goal);
	while (true) {
		ll from = from_list[goal];
		path.emplace_back(from);
		if (from == start) {
			break;
		}
		goal = from;
	}
	reverse(all(path));
	return path;
}



class FordFulkerson {
private:
	vb usedNode;
public:
	struct RevEdge { ll from, to, cap, rev; };

	FordFulkerson(ll n, Graph graph) 
		:usedNode(vb(n)), G(vec_t<2,RevEdge>(n))
	{
		rep(i, 0, graph.size()) {
			for (ll eind : graph.out_edges[i]) {
				Edge& e = graph.edges[eind];
				add_revedge(e);
			}
		}

	}
	vec_t<2, RevEdge> G;
	void add_revedge(Edge e) {
		G[e.from].push_back(RevEdge{ e.from, e.to ,e.cost, SZ(G[e.to]) });
		G[e.to].push_back(RevEdge{ e.to, e.from, 0 , SZ(G[e.from]) - 1 });
	}

	ll single_flow(ll from, ll to, ll flow) {
		// fromからtoに向かってflowを超えない範囲で一本のFlowを流す。
		if (from == to)
			return flow;
		usedNode[from] = 1;
		rep(i, 0, G[from].size()) {
			RevEdge& e = G[from][i];
			if (usedNode[e.to] || e.cap <= 0)
				continue;
			ll flow_from_e = single_flow(e.to, to, min(flow, e.cap));
			if (flow_from_e > 0) {
				e.cap -= flow_from_e; assert(e.cap >= 0);
				G[e.to][e.rev].cap += flow_from_e;
				// 今までよりも最大流を増やすことに成功したのでrerurn
				return flow_from_e;
			}
		}
		//すでにfromから先すべての辺を訪れていたあるいはすべてのcapが0だったら流せない。
		return 0;
	}
	ll max_flow(ll from, ll to) {
		ll flow = 0;
		while (true) {
			fill_v(usedNode, 0);
			ll f = single_flow(from, to, INF);
			if (f == 0)
				return flow;
			else
				flow += f;
		}

	}

};


// ================= Rectangle Area Problem =====================
auto getNeighbor = [](ll i, ll w, ll h) {
	ll H = i / w;
	ll W = i % w;
	vll res;
	if (H > 0) res.push_back(i - w);
	if (H < h - 1) res.push_back(i + w);
	if (W > 0)res.push_back(i - 1);
	if (W < w - 1)res.push_back(i + 1);
	return res;
};

auto getHW = [](ll i, ll w) {
	ll H = i / w;
	ll W = i % w;
	return pll{ H,W };
};













// ============================ Header  =================================

int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll N;
	cin >> N;
	vector<double> p(N);
	rep(i, 0, N)cin >> p[i];
	auto dp = make_v<double>(N + 1, N + 1);
	dp[0][0] = 1;
	rep(i, 0, N)rep(j, 0, N ) {
		dp[i + 1][j + 1] += dp[i][j]*p[j];
		dp[i][j+1] += dp[i][j] * (1-p[j]);
	}
	double res = 0;
	rep(i, 0, N + 1) {
		if (i > N / 2)
			res += dp[i][N];
	}
	cout << res << endl;
	return 0;

}

