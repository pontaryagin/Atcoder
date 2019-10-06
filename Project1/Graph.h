#pragma once
#include "MyHeader.h"
#include "UnionFind.h"
#include <fstream>
#ifdef _WIN64
	#include <direct.h>
#endif


enum class GraphDir { dir, undir };

template<class cost_t>
struct Edge_Base
{
	ll from;
	ll to;
	cost_t cost;
	Edge_Base reverse() const { return Edge_Base{ to, from , cost }; }
	Edge_Base(ll from , ll to, cost_t cost=1) : from(from),to(to),cost(cost){};
	Edge_Base(pll e) :from(e.first), to(e.second), cost(1) { }
	Edge_Base() :from(0), to(0), cost(1){ };
	bool operator<  (const Edge_Base& e) const {	return cost < e.cost; }
	bool operator>  (const Edge_Base& e) const {	return cost > e.cost; }
	bool operator== (const Edge_Base& e) const { return cost == e.cost && from == e.from && to == e.to; }
	string dot(GraphDir dir, bool weighted) const {
		return to_string(from) + (dir == GraphDir::dir ? "->": "--") + to_string(to) +
			(weighted? "[label = " + to_string(cost) + "]" :"" );
	}
};
using Edge = Edge_Base<ll>;

template<class EdgeType, class EdgeContainerType>
struct Edge_Itr_Base {
	constexpr Edge_Itr_Base() :index(), edges(nullptr) {}
	constexpr Edge_Itr_Base(ll index, EdgeContainerType& edges_) :index(index), edges(&edges_) {}
	constexpr Edge_Itr_Base& operator++() { ++index; return *this; }
	constexpr bool operator==(const Edge_Itr_Base& rhs) const { return index == rhs.index; }
	constexpr bool operator!=(const Edge_Itr_Base& rhs) const { return index != rhs.index; }
	constexpr EdgeType* operator->() const { return &(*edges)[index]; }
	constexpr EdgeType& operator*() const { return (*edges)[index]; }
	constexpr Edge_Itr_Base& operator+=(ll n) { index += n; return *this; }
	ll index;
	EdgeContainerType* edges;
};

auto nullAction = [](const auto&) {};

template<GraphDir dir, class cost_t>
struct Graph_Base {
	using Edge = Edge_Base<cost_t>;
	using Edge_Itr = Edge_Itr_Base<Edge_Base<cost_t>, vector<Edge_Base<cost_t>>>;
	using Edge_CItr = Edge_Itr_Base<const Edge_Base<cost_t>, const vector<Edge_Base<cost_t>>>;
	ll nodeSize;
	vector<Edge> edges;
	vector<vector<Edge_Itr>> out_edges;
	vector<vector<Edge_Itr>> in_edges;
	Graph_Base(ll nodeSize, const vector<Edge>& edges_ = vector<Edge>(), GraphDir dirct= GraphDir::dir)
		: nodeSize(nodeSize), out_edges(nodeSize), in_edges(nodeSize){
		for (const Edge& e : edges_) push(e);
	}
	Graph_Base(ll nodeSize, vector<pll> edges_)
		: nodeSize(nodeSize), out_edges(nodeSize), in_edges(nodeSize){
		for (const pll& e : edges_) push(Edge(e));
	}
	Graph_Base(vvll ajacency_matrix, ll default_value) 
		: nodeSize(ajacency_matrix.size()), out_edges(nodeSize), in_edges(nodeSize){
		ll n = ajacency_matrix.size();
		rep(i, 0, n)rep(j, 0, n) {
			if (ajacency_matrix[i][j] != default_value)
				push(Edge(i, j, ajacency_matrix[i][j]));
		}
	}
	Graph_Base(const Graph_Base& g) : nodeSize(g.nodeSize), out_edges(nodeSize), in_edges(nodeSize) {
		this->push(g.edges);
	}
	Graph_Base& operator=(const Graph_Base& g) {
		nodeSize = g.nodeSize; out_edges.resize(nodeSize); in_edges.resize(nodeSize);
		push(g.edges);
		return *this;
	}

	Edge& operator[](ll ind) { return this->edges[ind]; } 
	const Edge& operator[](ll ind) const{ return this->edges[ind]; }
	vector<Edge_Itr>& out(ll ind){ return this->out_edges[ind]; }
	const vector<Edge_Itr>& out(ll ind) const { return this->out_edges[ind]; }
	vector<Edge_Itr>& in(ll ind){ return this->in_edges[ind]; }
	const vector<Edge_Itr>& in(ll ind) const{ return this->in_edges[ind]; }
	Edge_Itr begin() { return Edge_Itr(0, edges); }
	Edge_Itr end() { return Edge_Itr(edges.size(), edges); }
	Edge_CItr begin() const { return Edge_CItr(0, edges); }
	Edge_CItr end() const { return Edge_CItr(edges.size(), edges); }

	ll size() const { return nodeSize; }
	ll sizeEdges() const { return edges.size(); }
	void _push(const Edge& edge){
		assert(max(edge.from, edge.to) < nodeSize);
		edges.emplace_back(edge);
		out_edges[edge.from].emplace_back(Edge_Itr(edges.size()-1,edges));
		in_edges[edge.to].emplace_back(Edge_Itr(edges.size() - 1, edges));
	}
public:
	template<class T = void>
	void push(const Edge& edge, enable_if_t<dir == GraphDir::undir, T*> = nullptr) {
		_push(edge); _push(edge.reverse());
	}
	template<class T = void>
	void push(const Edge& edge, enable_if_t<dir == GraphDir::dir, T*> = nullptr) {
		_push(edge);
	}
	void push(vector<Edge> edges) {
		for (const Edge& e : edges) {
			push(e);
		}
	}
	vvll adjacency_matrix() const {
		vvll d(size(), vll(size(), INF));
		for (auto& e : edges) {
			d[e.from][e.to] = e.cost;
		}
		rep(i, 0, size()) {
			d[i][i] = 0;
		}
		return d;
	}

	Graph_Base reverse() const {
		Graph_Base g(size());
		for (const auto& e : this->edges) {
			g.push(e.reverse());
		}
		return g;
	}

	vll get_topologically_sorted_nodes() const
	{
		// graph needs to be represented by adjacent list.
		// complexity: O( node size + edge size)
		ll nodeSize = this->size();

		// find root
		vll roots;
		vll inDegree(nodeSize);
		rep(i, 0, nodeSize)
		{
			for (auto& sibling : this->out(i)) {
				inDegree[sibling->to]++;
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
			for (auto& sibling : this->out(parent)) {
				inDegree[sibling->to]--;
				if (inDegree[sibling->to] == 0) {
					parents.push(sibling->to);
				}
			}
		}
		return sortedNodes;
	}

	// not safe for Edge_Itr
	//void topological_sort() {
	//	vll sorted = get_topologically_sorted_nodes();
	//	vll new_ind(sorted.size());
	//	vector<Edge> new_edges;
	//	rep(i, 0, sorted.size()) {
	//		new_ind[sorted[i]] = i;
	//	}
	//	for (Edge& e : edges) {
	//		new_edges.emplace_back(Edge{ new_ind[e.from], new_ind[e.to],e.cost });
	//	}
	//	*this = Graph_Base(this->size(), new_edges);
	//}
	cost_t diameter() const
	{
		// require : graph is tree
		// calculate the diameter ( longest path length ) in O(N)
		vector<cost_t> dp(size(), -1);
		cost_t m = 0; ll ind;
		function<void(ll)> dfs = [&](ll x) {
			for (auto& e : this->out(x)) {
				ll nextnode = e->to;
				if (dp[nextnode] == -1) {
					dp[nextnode] = dp[x] + e->cost;
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

	cost_t max_length() const {
		// calculate the max lenth of path in the graph
		vector<cost_t> dp(size(), -1);
		auto dfs = [&](auto dfs, ll x) -> void{
			for (auto& e : this->out(x)) {
				if (dp[e->to] == -1) {
					dp[e->to] = 0;
					dfs(dfs, e->to);
				}
				chmax(dp[e->from], dp[e->to] + e->cost);
			}
		};
		rep(node, 0, size())
			dfs(dfs, node);
		return *max_element(all(dp));
	}

	vll leaves() const {
		vll res;
		rep(i, 0, nodeSize) {
			if (out(i).size() <= 1)
				res.push_back(i);
		}
		return res;
	}

	template<class T, class S = decltype(nullAction)>
	void dfs(ll startNode, T before_act, S after_act = nullAction) const
	{
		// Impliment func: void(const Edge&) representing what this should do, when target node moves from visited node (e.from) to unvisited node (e.to).
		const auto& graph = *this;
		vb visited(graph.size());
		auto dfs_impl = [&](auto dfs_impl, ll startNode)-> void {
			visited[startNode] = 1;
			for (auto& e : graph.out(startNode)) {
				if (visited[e->to])
					continue;
				before_act(*e);
				dfs_impl(dfs_impl, e->to);
				after_act(*e);
			}
		};
		dfs_impl(dfs_impl, startNode);

	};

	template<class T, class S = decltype(nullAction)>
	void dfs_node(ll startNode, T before_act, S after_act = nullAction) const
	{
		// Impliment func: void(ll node_ind) representing what this should do, when target node moves from visited node to unvisited node (node_ind).
		const auto& graph = *this;
		vb visited(graph.size());
		auto dfs_impl = [&](auto dfs_impl, ll startNode)-> void {
			before_act(startNode);
			visited[startNode] = 1;
			for (auto& e : graph.out(startNode)) {
				if (visited[e->to])
					continue;
				dfs_impl(dfs_impl, e->to);
			}
			after_act(startNode);
		};
		dfs_impl(dfs_impl, startNode);
	};

	template<class T, class S = decltype(nullAction)>
	void bfs(ll startNode, T before_act, S after_act = nullAction) const
	{
		const auto& graph = *this;
		vb visited(graph.size());
		auto bfs_impl = [&](ll startNode) {
			//if (visited[startNode] != 0) return;
			visited[startNode] = 1;
			queue<Edge> toVisit;
			for (auto& e : graph.out(startNode))
				toVisit.push(*e);
			while (toVisit.size()) {
				auto next = toVisit.front(); toVisit.pop();
				if (visited[next.to])
					continue;
				visited[next.to] = 1;
				before_act(next);
				for (auto& e : graph.out(next.to)) {
					if (!visited[e->to])
						toVisit.push(*e);
				}
				after_act(next);
			}
		};
		bfs_impl(startNode);
	};
	
	vector<cost_t> dijkstra(ll start) const {
		vector<cost_t> fromList;
		return dijkstra(start, fromList);
	}

	vector<cost_t> dijkstra(ll start, vector<cost_t>& from_list) const {
		// graph: weighted directed graph of adjacent representation
		// start: index of start point
		// return1: minimum path length from start
		// complexity : E*log(V)
		const auto& graph = *this;
		ll node_size = graph.size();
		vector<cost_t> dist(node_size, INF);
		from_list.resize(node_size);
		fill_v(from_list, -1);
		dist[start] = 0;
		pq_greater<pair<cost_t, pll>> pq;
		pq.push({ 0, {start, start} });
		while (!pq.empty()) {
			auto node = pq.top(); pq.pop();
			// if not shortest path fixed, fix
			ll from = node.second.first;
			ll to = node.second.second;
			if (from_list[to] != -1)
				continue;
			from_list[to] = from;

			for (auto& edge : graph.out(to)) {
				ll adj = edge->to;
				cost_t cost = dist[to] + edge->cost;
				if (dist[adj] > cost) {
					dist[adj] = min(dist[adj], cost);
					pq.push({ cost ,{to, adj} });
				}
			}
		}
		return dist;
	}

	vll euler_tour(ll start) const
	{
		vll res;
		res.push_back(start);
		dfs(start, [&](const Edge& e) {
			res.push_back(e.to);
			}, [&](const Edge& e) {
				res.push_back(e.from);
			});
		return res;
	}

	template<GraphDir out_dir>
	Graph_Base<out_dir, cost_t> kruskal() const
	{
		//returns minimal spanning tree
		Graph_Base<out_dir, cost_t> res(nodeSize);
		vpll sortedEdges;
		rep(i, 0, edges.size()) {
			sortedEdges.push_back({ edges[i].cost, i });
		}
		sort(all(sortedEdges));
		UnionFind uf(nodeSize);
		rep(i, 0, sortedEdges.size()) {
			ll cost, eInd;
			tie(cost, eInd) = sortedEdges[i];
			ll from = (*this)[eInd].from; ll to = (*this)[eInd].to;
			if (!uf.is_same(from, to)) {
				res.push((*this)[eInd]);
			}
			uf.unite(from, to);
		}
		return res;
	}

	vvll warshall_floyd() const {
		// O(|V|^3)
		const Graph_Base& g = *this;
		ll n = g.size();
		vvll d = g.adjacency_matrix();
		rep(k, 0, n)rep(i, 0, n)rep(j, 0, n) {
			if (d[i][j] > d[i][k] + d[k][j])
				d[i][j] = d[i][k] + d[k][j];
		}
		return d;
	}

	vll bellman_ford(ll start, ll negative_closed_loop_value = -INF) const {
		vll from_list;
		return bellman_ford(start, from_list, negative_closed_loop_value);
	}

	vll bellman_ford(ll start, vll& from_list, ll negative_closed_loop_value = -INF) const {
		// O(|E| * |V|)
		const Graph_Base& g = *this;
		vll dist(g.size(), INF);
		dist[start] = 0;
		from_list.resize(g.size());
		rep(i, 0, g.size()) {
			for(const Edge& e: g){
				if (dist[e.from] != INF && dist[e.to] > dist[e.from] + e.cost) {
					dist[e.to] = dist[e.from] + e.cost;
					from_list[e.to] = e.from;
					if (i == g.size() - 1 && dist[e.to] != INF) {
						// check negative closed loop
						dist[e.to] = negative_closed_loop_value;
					}
				}
			}
		}
		// propagate negative path
		rep(i, 0, g.size()) {
			rep(j, 0, g.edges.size()) {
				auto& e = g.edges[j];
				if (dist[e.from] == negative_closed_loop_value && dist[e.from] != INF)
					dist[e.to] = negative_closed_loop_value;
			}
		}
		return dist;
	}

	bool is_bipartite() const {
		vll even(size(),-1);
		even[0] = 0;
		bool ok = true;
		dfs_node(0,
			[&](ll node) {
				for (auto& e : out(node)) {
					if (even[e->to] != -1 ) {
						if (even[e->from] == even[e->to]) {
							ok = false;
							break;
						}
					}
					else {
						even[e->to] = !even[e->from];
					}
				}
			});
		return ok;
	}


	bool acyclic() const {
		vll loop; 
		return acyclic(loop);
	}
	bool acyclic(vll& loop) const {
		// check whether directed graph has cycle in O(|V| + |E|)
		// found loop is stored to "loop"
		auto& g = *this;
		vll visited(size());
		vll stack;
		auto dfs = [&](auto dfs, ll node, ll par) -> bool {
			visited[node] = 1; stack.push_back(node);
			for (auto& e : g.out(node)) {
				if ((dir == GraphDir::dir || par != e->to) && visited[e->to] == 1) {
					if (loop.empty()) {
						auto it = find(all(stack), e->to);
						copy(it, stack.end(), back_inserter(loop));
					}
					return false;
				}
				else if (visited[e->to] == 0 && !dfs(dfs, e->to, e->from))
					return false;
			}
			visited[node] = 2; stack.pop_back();
			return true;
		};
		rep(i, 0, g.size()) {
			if (!dfs(dfs, i, i))
				return false;
		}
		return true;
	}

	Graph_Base scc(vll& components) const {
		// strongly connected components decomposition algorithm in O(|V| + |E|)
		// @ return : contracted DAG
		// @ in (components) : node i is in components[i]-th component in DAG
		static_assert(dir == GraphDir::dir, "scc is valid for directed graph");
		components.resize(size());
		vpll time(size());
		ll now = 0;
		vll visited(size());
		auto dfs = [&](auto dfs, ll node) -> void {
			visited[node] = 1; 
			for (auto&& e : this->out(node)) {
				if (visited[e->to] == 0) {
					dfs(dfs, e->to);
				}
			}
			time[node] = { now++, node };
		};
		rep(i, 0, size()) {
			if (!visited[i]) {
				dfs(dfs, i);
			}
		}
		fill_v(visited, 0);
		sort(all(time));
		ll cur_comp = 0;
		auto dfs_rev = [&](auto dfs_rev, ll node) -> void {
			visited[node] = 1; components[node] = cur_comp;
			for (auto&& e : this->in(node)) {
				if (visited[e->from] == 0) {
					dfs_rev(dfs_rev, e->from);
				}
			}
		};
		rrep(i, 0, time.size()) {
			ll node = time[i].second;
			if(visited[node]== 0)
				dfs_rev(dfs_rev, node);
			cur_comp++;

		}
		// create contracted DAG 
		Graph_Base res(cur_comp);
		set<pll> check;
		for(const Edge& e: edges) {
			ll from = components[e.from]; ll to = components[e.to];
			if (from != to && exist(check, {from, to})) {
				res.push({ from, to});
			}
		}
		return res;

	}
	string dot(bool weighted = false) const {
		// export graph as dot file
		string res = (dir == GraphDir::dir?"digraph {" : "graph {");
		res += "node [ style=filled shape=circle fontname=\"Fira Code\" fontcolor=darkslategray color=darkslategray fillcolor=lightcyan];\n";
		res += "edge [ color=darkslategray ]";
		rep(i, 0, size()) {
			res += to_string(i) + ";\n";
		}
		set<pll> used; 
		for (auto& e : edges) {
			if (!exist(used, { e.from, e.to }) && (dir == GraphDir::dir || !exist(used, { e.to, e.from }))) {
				used.insert({ e.from, e.to});
				res += e.dot(dir, weighted);
				res += ";\n";
			}
		}
		res += "}\n";
		return res;
	}
	void show(bool weighted = false, string  ext = "pdf") const {
		// show graph as png file
#ifdef _WIN64
		srand(time(nullptr));
		(void)_mkdir("./tmp");
		string tmpdot = "./tmp/"; string tmppng = "./tmp/";
		tmpdot += to_string(rand()); tmppng += to_string(rand()) + "." + ext;
		ofstream of(tmpdot);
		of <<  dot(weighted) << endl;
		(int)system(("dot -T"+ext+ " -o "+ tmppng + " " + tmpdot).c_str());
		(int)system(("powershell start " + tmppng).c_str());
#endif // _WIN64
	}
};

using Graph = Graph_Base<GraphDir::undir, ll>;
using Digraph = Graph_Base<GraphDir::dir, ll>;


vll shortest_path_generator(const vll& from_list, ll start, ll goal) {
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
	using Edge = Edge_Base<ll>;
public:
	struct RevEdge { ll from, to, cap, rev; };
	FordFulkerson(Digraph graph) 
		:usedNode(graph.size()), G(vec_t<2,RevEdge>(graph.size()))
	{
		rep(i, 0, graph.size()) {
			for (auto& e : graph.out(i)) {
				add_revedge(*e);
			}
		}

	}
	vec_t<2, RevEdge> G;
	void add_revedge(const Edge& e) {
		G[e.from].push_back(RevEdge{ e.from, e.to ,e.cost, SZ(G[e.to]) });
		G[e.to].push_back(RevEdge{ e.to, e.from, 0 , SZ(G[e.from]) - 1 });
	}
	
	ll single_flow(ll from, ll to, ll flow) {
		// make a single flow
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
				// get a larger flow
				return flow_from_e;
			}
		}
		// if we already visited all edges or cap = 0 flow = 0;
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

// Least Common Ancestor
class LCA {
public:
	LCA(const Graph& graph, ll root) : max_par(ll(ceil(log2(graph.size()) + 2))), parent(graph.size(), vll(max_par,-1)),
		depth() {
		//parent[root][0] = root;
		graph.dfs(root, [&](const Edge & e) {
			ll to = e.to;
			parent[to][0] = e.from;
			rep(i, 1, parent[to].size()) {
				if (parent[to][i - 1] == -1)
					return;
				else
					parent[to][i] = (parent[parent[to][i - 1]][i - 1]);
			}
		});
		depth = graph.dijkstra(root);
	}
	ll operator()(ll node1, ll node2) {
		if (depth[node1] > depth[node2]) swap(node1, node2);
		rrep(i, 0, max_par) {
			if (((depth[node2] - depth[node1]) >> i) & 1) {
				node2 = parent[node2][i];
			}
		}
		if (node1 == node2)return node1;
		rrep(i, 0, max_par) {
			if (parent[node1][i] != parent[node2][i]) {
				node1 = parent[node1][i];
				node2 = parent[node2][i];
			}
		}
		return parent[node1][0];
	}
private:
	ll max_par;
	vvll parent;
	vll depth;
};

class BipartiteMatching {
	// O(V*E)
	int n, left, right;
	vector< vector< int > > graph;
	vector< int > used;
	int timestamp;
public:
	BipartiteMatching(int left, int right) : n(left+right), left(left), right(right), graph(n), used(n, 0), timestamp(0){}

	void push(int u, int v) {
		graph[u].push_back(v + left);
		graph[size_t(v) + left].push_back(u);
	}

	bool dfs(int idx, vector<int>& match) {
		used[idx] = timestamp;
		for (auto& to : graph[idx]) {
			int to_match = match[to];
			if (to_match == -1 || (used[to_match] != timestamp && dfs(to_match, match))) {
				match[idx] = to;
				match[to] = idx;
				return true;
			}
		}
		return false;
	}

	int bipartite_match(vector<int>& match) {
		match.resize(n); fill_v(match, -1);
		int ret = 0;
		for (int i = 0; i < SZ(graph); i++) {
			if (match[i] == -1) {
				++timestamp;
				ret += dfs(i, match);
			}
		}
		return ret;
	}

	int bipartite_match() {
		vector<int> match;
		return bipartite_match(match);
	}
};

// tree dfs
//auto dfs = [&](auto dfs, ll from, ll to) ->void {
//	for (auto& e : g.out(to)) {
//		if (e->to != from) {
//			// do something
//		}
//	}
//};
//dfs(dfs, 0, 0);