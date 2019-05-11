#pragma once
#include "MyHeader.h"

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
	enum Dir{dir, nondir};
	Graph(ll nodeSize, const vector<Edge>& edges = vector<Edge>(), Dir dirct= dir)
		: nodeSize(nodeSize), out_edges(nodeSize), in_edges(nodeSize), edges(move(edges)){
		rep(i, 0, edges.size()) {
			in_edges[edges[i].to].push_back(i);
			out_edges[edges[i].from].push_back(i);
			if (dirct == nondir) {
				in_edges[edges[i].from].push_back(i);
				out_edges[edges[i].to].push_back(i);
			}
		}
	}
	Graph(ll nodeSize, vector<pll> edges_, Dir dirct = dir)
		: nodeSize(nodeSize), out_edges(nodeSize), in_edges(nodeSize), edges(edges_.size()) {
		//if (dirct == nondir) edges.resize(edges_.size() * 2);
		rep(i, 0, edges.size()) {
			edges[i] = edges_[i];
			in_edges[edges[i].to].push_back(i);
			out_edges[edges[i].from].push_back(i);
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
	
	Edge& operator[](ll ind) {
		return this->edges[ind];
	}
	const Edge& operator[](ll ind) const {
		return this->edges[ind];
	}
	size_t size() const { return out_edges.size(); }
	void push(Edge edge){
		assert(max(edge.from, edge.to) < out_edges.size());
		edges.emplace_back(edge);
		out_edges[edge.from].emplace_back(edges.size() - 1);
		in_edges[edge.to].emplace_back(edges.size() - 1);
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
	ll diameter() {
		// require : graph is tree
		// calculate the diameter ( longest path length ) in O(N)
		vll dp(size(),-1);
		ll m = 0; ll ind;
		function<void(ll)> dfs = [&](ll x) {
			for (ll e : out_edges[x]) {
				ll nextnode = edges[e].to;
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
