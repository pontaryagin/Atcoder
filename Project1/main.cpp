
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

template<typename Inputs>
void sort_uniq(Inputs& inputs) {
	sort(all(inputs));
	inputs.erase(unique(all(inputs)), inputs.end());
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
	bool operator<  (const Edge& e) const {
		return cost < e.cost;
	}
	bool operator>  (const Edge& e) const {
		return cost > e.cost;
	}
};


struct Graph {
	ll nodeSize;
	vector<Edge> edges;
	vector<vector<ll>> out_edges;
	vector<vector<ll>> in_edges;
	enum Dir{dir, undir};
	Graph(ll nodeSize, const vector<Edge>& edges_ = vector<Edge>(), Dir dirct= dir)
		: nodeSize(nodeSize), out_edges(nodeSize), in_edges(nodeSize), edges(){
		if (dirct == undir) {
			for (const Edge& e : edges_) push_undir(e);
		}
		else {
			for (const Edge& e : edges_) push(e);
		}

	}
	Graph(ll nodeSize, vector<pll> edges_, Dir dirct = dir)
		: nodeSize(nodeSize), out_edges(nodeSize), in_edges(nodeSize), edges() {
		if (dirct == undir) {
			for (const pll& e : edges_) push_undir(Edge(e));
		}
		else {
			for (const pll& e : edges_) push(Edge(e));
		}
	}
	Graph(vvll ajacency_matrix, ll default_value) 
		: nodeSize(ajacency_matrix.size()), out_edges(nodeSize), in_edges(nodeSize){
		ll n = ajacency_matrix.size();
		rep(i, 0, n)rep(j, 0, n) {
			if (ajacency_matrix[i][j] != default_value)
				push(Edge(i, j, ajacency_matrix[i][j]));
		}
	}
	
	Edge& operator[](ll ind) { return this->edges[ind]; } 
	const Edge& operator[](ll ind) const{ return this->edges[ind]; }
	vector<ll>& out(ll ind){ return this->out_edges[ind]; }
	const vector<ll>& out(ll ind) const { return this->out_edges[ind]; }
	vector<ll>& in(ll ind){ return this->in_edges[ind]; }
	const vector<ll>& in(ll ind) const{ return this->in_edges[ind]; }

	size_t size() const { return nodeSize; }

	void push(const Edge& edge){
		assert(max(edge.from, edge.to) < nodeSize);
		edges.emplace_back(edge);
		out_edges[edge.from].emplace_back(edges.size()-1);
		in_edges[edge.to].emplace_back(edges.size() - 1);
	}
	void push_undir(const Edge& edge) {
		push(edge); push(edge.reverse());
	}

	//void erase(ll from, ll to) {
	//	// O(nodeSize)
	//	assert(max(from, to) < out_edges.size());
	//	for (Edge& e : out_edges[from]) {
	//		if (e.to = to) {
	//			edges.erase(edges.begin() + k);
	//			out_edges
	//		}
	//	}
	//}
	//void erase(ll edge_ind) {
	//	// O(nodeSize)
	//	assert(edge_ind < edges.size());
	//	edges.erase(edges.begin() + edge_ind);
	//	auto& esin = in_edges[edges[edge_ind].to];
	//	esin.erase(find(all(esin), edge_ind));
	//	auto& esout = out_edges[edges[edge_ind].from];
	//	esout.erase(find(all(esout), edge_ind));
	//}
	void push(vector<Edge> edges) {
		for (const Edge& e : edges) {
			push(e);
		}
	}
	vvll adjacency_matrix(ll default_value = INF) const {
		vvll d(size(), vll(size()));
		for (auto& e : edges) {
			d[e.from][e.to] = e.cost;
		}
		return d;
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
			for (ll sibling_ind : this->out(i)) {
				inDegree[(*this)[sibling_ind].to]++;
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
			for (ll sibling_ind : this->out(parent)) {
				auto& sibling = (*this)[sibling_ind];
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
	ll diameter() {
		// require : graph is tree
		// calculate the diameter ( longest path length ) in O(N)
		vll dp(size(),-1);
		ll m = 0; ll ind;
		function<void(ll)> dfs = [&](ll x) {
			for (ll e_ind : out(x)) {
				ll nextnode = (*this)[e_ind].to;
				if (dp[nextnode] == -1) {
					dp[nextnode] = dp[x] + 1;
					if (dp[nextnode] > m) {
						m = dp[nextnode];  ind = nextnode;
					}
					dfs(nextnode);
				}
			}
		};
		dp[0] = 0; ind = 0;
		dfs(0);
		ll first = ind;
		fill_v(dp, -1);
		dp[first] = 0; 
		dfs(first);
		return m;
		// remark two end points of diameter are 'first' and 'ind';
	}

	vll leaves() {
		vll res;
		rep(i, 0, nodeSize) {
			if (out(i).size() <= 1)
				res.push_back(i);
		}
		return res;
	}

};

class GraphDFS
{
	vb visited;
	Graph* graph;
public:
	GraphDFS(Graph& graph) :visited(graph.size()), graph(&(graph)){}

	// Impliment func: void(Edge&) representing what this should do, when target node moves from visited node (e.from) to unvisited node (e.to).
	template<class T>
	void operator()(T&& func, ll startNode=0) {

		if (visited[startNode] != 0) return;
		visited[startNode] = 1;
		for (ll e_ind: graph->out(startNode)) {
			auto& e = (*graph)[e_ind];
			func(e);
			operator()(func, e.to);
		}
	}
};

pair<vll, vll> dijkstra(const Graph& graph, size_t start) {
	// graph: weighted directed graph of adjacent representation
	// start: index of start point
	// return1: minimum path length from start
	// return2: concrete shortest path info
	// complexity : E*log(V)
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

		for (ll edge_ind : graph.out(to)) {
			auto& edge = graph[edge_ind];
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

pair<vll, vll> dijkstra2(const Graph& graph, size_t start) {
	// graph: weighted directed graph of adjacent representation
	// start: index of start point
	// return1: minimum path length from start
	// return2: concrete shortest path info
	// complexity : E*log(V)
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

		for (ll edge_ind : graph.out(to)) {
			auto& edge = graph[edge_ind];
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

vvll warshall_floyd(const Graph& g, ll default_value) {
	ll n = g.size();
	vvll d = g.adjacency_matrix(INF);
	rep(k, 0, n)rep(i, 0, n)rep(j, 0, n) {
		if (d[i][j] > d[i][k] + d[k][j])
			d[i][j] = d[i][k] + d[k][j];
	}
	rep(i, 0, n)rep(j, 0, n) {
		if (d[i][j] == INF)
			d[i][j] = default_value;
	}
	return d;
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
			for (ll e_ind : graph.out(i)) {
				add_revedge(graph[e_ind]);
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




struct UnionFind {
	vector<ll> data;
	vll querySize_;
	set<ll> roots;
	UnionFind(ll size) : data(size, -1), querySize_(size, 0) {
		rep(i, 0, size) roots.insert(i);
	}

	ll unite(ll x, ll y) {
		// return: root
		x = operator[](x); y = operator[](y);
		if (x != y) {
			if (data[y] < data[x]) swap(x, y);
			data[x] += data[y]; data[y] = x;
			querySize_[x] += querySize_[y] + 1;
			roots.erase(y);
			return x;
		}
		else {
			querySize_[x]++;
			return x;
		}
	}
	bool is_same(ll x, ll y) {
		// check whether x and y are connected
		return operator[](x) == operator[](y);
	}
	ll operator[](ll x) {
		// get root
		return data[x] < 0 ? x : data[x] = operator[](data[x]);
	}
	ll size(ll x) {
		return -data[operator[](x)];
	}
	ll  query_size(ll x) {
		return querySize_[operator[](x)];
	}
	const set<ll>& getRoots() {
		return roots;
	}
	ll rank(ll x) {
		return -data[operator[](x)];
	}
	void initialize() {
		for (auto& i : data) {
			i = -1;
		}
	}
};












// ============================ Header  =================================



int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	ll n, m; cin >> n >> m;
	vpll ab(m);
	read(ab);

	UnionFind uf(n);
	Graph g(n);
	rep(i, 0, m) {
		uf.unite(ab[i].fir-1, ab[i].sec-1);
		g.push(Edge{ ab[i].first - 1,ab[i].sec - 1 });

	}
	bool ok = 1;
	rep(i, 0, n) {
		if ((g.out(i).size() + g.in(i).size()) % 2) {
			ok = 0;
		}
	}
	if (ok)
		cout << "YES" << endl;
	else
		cout << "NO" << endl;

	return 0;

}
