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