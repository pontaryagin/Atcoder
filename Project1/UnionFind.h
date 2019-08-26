#pragma once
#include "MyHeader.h"

struct UnionFind {
	vector<ll> data;
	vll querySize_;
	set<ll> roots;
	UnionFind(ll size) : data(size, -1), querySize_(size, 0) {
		rep(i, 0, size) roots.insert(i);
	}

	ll unite(ll x, ll y) {
		// return: root
		x = get_root(x); y = get_root(y);
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
		return get_root(x) == get_root(y);
	}
	ll get_root(ll x) {
		// get root
		return data[x] < 0 ? x : data[x] = get_root(data[x]);
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
};