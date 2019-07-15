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

auto nullAction = [](const Edge&) {};

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

	vll get_topologically_sorted_nodes() const;

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
	ll diameter() const;

	vll leaves() const {
		vll res;
		rep(i, 0, nodeSize) {
			if (out(i).size() <= 1)
				res.push_back(i);
		}
		return res;
	}
	template<class T, class S = decltype(nullAction)>
	void dfs(ll startNode, T before_act, S after_act = nullAction) const;
	// Impliment func: void(Edge&) representing what this should do, when target node moves from visited node (e.from) to unvisited node (e.to).
	template<class T, class S = decltype(nullAction)>
	void bfs(ll startNode, T before_act, S after_act = nullAction) const;
	// Impliment func: void(Edge&) representing what this should do, when target node moves from visited node (e.from) to unvisited node (e.to).
	vll dijkstra(ll start) const;
	vll dijkstra(ll start, vll& from_list) const;
	vll euler_tour(ll start) const;
};

vll Graph::get_topologically_sorted_nodes() const {
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

template<class T, class S>
void Graph::dfs(ll startNode, T beforeAct, S afterAct) const
{
	// Impliment func: void(const Edge&) representing what this should do, when target node moves from visited node (e.from) to unvisited node (e.to).
	const auto& graph = *this;
	vb visited(graph.size());
	auto dfs_impl = [&](auto dfs_impl, ll startNode)-> void {
		visited[startNode] = 1;
		for (ll e_ind : graph.out(startNode)) {
			auto& e = graph[e_ind];
			if (visited[e.to])
				continue;
			beforeAct(e);
			dfs_impl(dfs_impl, e.to);
			afterAct(e);
		}
	};
	dfs_impl(dfs_impl, startNode);

};

ll Graph::diameter() const {
	// require : graph is tree
	// calculate the diameter ( longest path length ) in O(N)
	vll dp(size(), -1);
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

template<class T, class S>
void Graph::bfs(ll startNode, T beforeAct, S afterAct) const
{
	const auto& graph = *this;
	vb visited(graph.size());
	auto bfs_impl = [&](auto bfs_impl, ll startNode) {
		//if (visited[startNode] != 0) return;
		visited[startNode] = 1;
		queue<Edge> toVisit;
		for (auto ei : graph.out(startNode))
			toVisit.push(graph[ei]);
		while (toVisit.size()) {
			auto next = toVisit.front(); toVisit.pop();
			if (visited[next.to])
				continue;
			visited[next.to] = 1;
			beforeAct(next);
			for (auto ei : graph.out(next.to)) {
				if (!visited[graph[ei].to])
					toVisit.push(graph[ei]);
			}
			afterAct(next);
		}
	};
	bfs_impl(bfs_impl, startNode);
};

vll Graph::dijkstra(ll start) const {
	vll fromList;
	return dijkstra(start, fromList);
}

vll Graph::dijkstra(ll start, vll& fromList) const {
	// graph: weighted directed graph of adjacent representation
	// start: index of start point
	// return1: minimum path length from start
	// complexity : E*log(V)
	const auto& graph = *this;
	ll node_size = graph.size();
	vll dist(node_size, 1LL << 60);
	fromList.resize(node_size);
	fill_v(fromList, -1);
	dist[start] = 0;
	pq_greater<pair<ll, pll>> pq;
	pq.push({ 0, {start, start} });
	while (!pq.empty()) {
		auto node = pq.top(); pq.pop();
		// if not shortest path fixed, fix
		ll from = node.second.first;
		ll to = node.second.second;
		if (fromList[to] != -1)
			continue;
		fromList[to] = from;

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
	return dist;
}

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

vll Graph::euler_tour(ll start) const {
	vll res;
	res.push_back(start);
	dfs(start, [&](const Edge& e) {
			res.push_back(e.to);
		}, [&](const Edge& e) {
			res.push_back(e.from);
		});
	//res.push_back(start);
	return res;
}


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
