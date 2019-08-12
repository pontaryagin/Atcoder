#pragma once
#include "MyHeader.h"
#include "UnionFind.h"

struct Edge
{
	ll from;
	ll to;
	ll cost=1;
	Edge reverse() const { return Edge{ to, from , cost }; }
	Edge(ll from , ll to, ll cost=1) : from(from),to(to),cost(cost){};
	Edge(pll e) { from = e.first; to = e.second; cost = 1; }
	Edge() :from(0), to(0), cost(0){ };
	bool operator<  (const Edge& e) const {	return cost < e.cost; }
	bool operator>  (const Edge& e) const {	return cost > e.cost; }
	bool operator== (const Edge & e) const { return cost == e.cost && from == e.from && to == e.to; }
};

struct Edge_Itr {
	Edge_Itr():index(), edges(nullptr){}
	Edge_Itr(ll index, vector<Edge>& edges_):index(index), edges(&edges_){}
	Edge_Itr& operator++() { ++index; return *this; }
	bool operator==(const Edge_Itr& rhs) const { return index == rhs.index; }
	Edge* operator->() const { return &(*edges)[index]; }
	Edge& operator*() const { return (*edges)[index]; }
	Edge_Itr& operator+=(ll n) { index += n; return *this; }
	ll index;
	vector<Edge>* edges;
};

auto nullAction = [](const Edge&) {};

struct Graph {
	ll nodeSize;
	vector<Edge> edges;
	vector<vector<Edge_Itr>> out_edges;
	vector<vector<Edge_Itr>> in_edges;
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
	vector<Edge_Itr>& out(ll ind){ return this->out_edges[ind]; }
	const vector<Edge_Itr>& out(ll ind) const { return this->out_edges[ind]; }
	vector<Edge_Itr>& in(ll ind){ return this->in_edges[ind]; }
	const vector<Edge_Itr>& in(ll ind) const{ return this->in_edges[ind]; }

	size_t size() const { return nodeSize; }

	void push(const Edge& edge){
		assert(max(edge.from, edge.to) < nodeSize);
		edges.emplace_back(edge);
		out_edges[edge.from].emplace_back(Edge_Itr(edges.size()-1,edges));
		in_edges[edge.to].emplace_back(Edge_Itr(edges.size() - 1, edges));
	}
	void push(const Edge& edge, Graph::Dir dir) {
		if (dir == Dir::undir)
			push_undir(edge);
		else
			push(edge);
	}
	void push_undir(const Edge& edge) {
		push(edge); push(edge.reverse());
	}
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
	ll diameter() const
	{
		// require : graph is tree
		// calculate the diameter ( longest path length ) in O(N)
		vll dp(size(), -1);
		ll m = 0; ll ind;
		function<void(ll)> dfs = [&](ll x) {
			for (auto& e : out(x)) {
				ll nextnode = e->to;
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
	void bfs(ll startNode, T before_act, S after_act = nullAction) const
	{
		const auto& graph = *this;
		vb visited(graph.size());
		auto bfs_impl = [&](auto bfs_impl, ll startNode) {
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
		bfs_impl(bfs_impl, startNode);
	};
	
	vll dijkstra(ll start) const {
		vll fromList;
		return dijkstra(start, fromList);
	}

	vll dijkstra(ll start, vll& from_list) const {
		// graph: weighted directed graph of adjacent representation
		// start: index of start point
		// return1: minimum path length from start
		// complexity : E*log(V)
		const auto& graph = *this;
		ll node_size = graph.size();
		vll dist(node_size, 1LL << 60);
		from_list.resize(node_size);
		fill_v(from_list, -1);
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

			for (auto& edge : graph.out(to)) {
				ll adj = edge->to;
				ll cost = dist[to] + edge->cost;
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

	Graph kruskal(Graph::Dir = Dir::undir) const
	{
		//returns minimal spanning tree
		Graph res(nodeSize);
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
				res.push((*this)[eInd], dir);
			}
			uf.unite(from, to);
		}
		return res;
	}

	vvll warshall_floyd(ll default_value = INF) const {
		// O(|V|^3)
		const Graph& g = *this;
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

	vll bellman_ford(ll start, ll negative_closed_loop_value = -INF) const {
		vll from_list;
		return bellman_ford(start, from_list, negative_closed_loop_value);
	}

	vll bellman_ford(ll start, vll& from_list, ll negative_closed_loop_value = -INF) const {
		// O(|E| * |V|)
		const Graph& g = *this;
		vll dist(g.size(), INF);
		dist[start] = 0;
		from_list.resize(g.size());
		rep(i, 0, g.size()) {
			rep(j, 0, g.edges.size()) {
				auto& e = g.edges[j];
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
};


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
public:
	struct RevEdge { ll from, to, cap, rev; };

	FordFulkerson(ll n, Graph graph) 
		:usedNode(vb(n)), G(vec_t<2,RevEdge>(n))
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


// Least Common Ancestor
class LCA {
public:
	LCA(const Graph& graph, ll root) : max_par(ceil(log2(graph.size()) + 2)), parent(graph.size(), vll(max_par,-1)),
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
