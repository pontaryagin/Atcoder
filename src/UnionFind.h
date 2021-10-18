#pragma once
#include "MyHeader.h"

struct UnionFind {
	vector<ll> data;
	vll querySize_;
	set<ll> roots;
	vll diff_weight;
	UnionFind(ll size) : data(size, -1), querySize_(size, 0), diff_weight(size, 0) {
		rep(i, 0, size) roots.insert(i);
	}
	// return : pair {new root, old root}
	pll unite(ll x, ll y, ll w = 0) {
		// return: root
		w += weight(x); w -= weight(y);
		x = get_root(x); y = get_root(y);
		if (x != y) {
			if (data[y] < data[x]) swap(x, y), w = -w;
			diff_weight[y] = w;
			data[x] += data[y]; data[y] = x;
			querySize_[x] += querySize_[y] + 1;
			roots.erase(y);
			return { x, y };
		}
		else {
			querySize_[x]++;
			return { x, y };
		}
	}
	bool is_same(ll x, ll y) {
		// check whether x and y are connected
		return get_root(x) == get_root(y);
	}
	ll get_root(ll x) {
		// get root and compress path
		if (data[x] < 0) {
			return x;
		} else {
			auto root = get_root(data[x]);
			diff_weight[x] += diff_weight[data[x]];
			return data[x] = root;
		}
	}
	ll size(ll x) {
		return -data[get_root(x)];
	}
	ll  query_size(ll x) {
		return querySize_[get_root(x)];
	}

	const set<ll>& get_roots() {
		return roots;
	}
	void initialize() {
		for (auto& i : data) {
			i = -1;
		}
	}
	ll weight(ll x) {
		get_root(x);
		return diff_weight[x];
	}
	ll weight_diff(ll x, ll y) {
		return weight(y) - weight(x);
	}
};