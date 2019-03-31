#pragma once
#include "MyHeader.h"

struct Edge
{
	ll from;
	ll to;
	ll cost;
};

using Graph = vector<vector<Edge>>;

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

		for (Edge edge : graph[to]) {
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

vll get_topologically_sorted_nodes(const vvll& graph)
{
	// graph needs to be represented by adjacent list.
	// complexity: O( node size + edge size)
	ll nodeSize = graph.size();

	// find root
	vll roots;
	vll inDegree(nodeSize);
	rep(i, 0, nodeSize)
	{
		for (ll sibling : graph[i]) {
			inDegree[sibling]++;
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
		for (ll sibling : graph[parent]) {
			inDegree[sibling]--;
			if (inDegree[sibling] == 0) {
				parents.push(sibling);
			}
		}
	}
	return sortedNodes;
}

// ================= Rectangle Area Problem =====================
auto getNeighbor = [&](ll i, ll w, ll h) {
	ll H = i / w;
	ll W = i % w;
	vll res;
	if (H > 0) res.push_back(i - w);
	if (H < h - 1) res.push_back(i + w);
	if (W > 0)res.push_back(i - 1);
	if (W < w - 1)res.push_back(i + 1);
	return res;
};

auto getHW = [&](ll i, ll w) {
	ll H = i / w;
	ll W = i % w;
	return pll{ H,W };
};
