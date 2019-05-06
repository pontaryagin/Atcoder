#include "MyHeader.h"
//#include "Algorithm.h"
//#include "Bit.h"
//#include "UnionFind.h"
//#include "NumberTheory.h"
#include "Graph.h"
//#include "SegmentTree.h"
//#include "Dijkstra.h"
//#include "bits/stdc++.h"
#include "boost/iterator/iterator_facade.hpp"

// ============================ Header  =================================

struct my_list_node {
	ll data;
	unique_ptr<my_list_node> next;
	my_list_node(ll data, unique_ptr<my_list_node>&& next) : data(data), next(move(next)) {	}
	my_list_node(my_list_node&& node) noexcept : data(node.data), next(move(node.next)) {	}
};



class my_lsit {
public:
	// define iter 
	class iterator
		:public boost::iterator_facade<iterator, ll, boost::forward_traversal_tag>
	{
		my_list_node* p_;
	public:
		iterator() = default;
		iterator(my_list_node* p) :p_(p) {}
	private:
		friend class boost::iterator_core_access;
		void increment() { p_ = p_->next.get(); }
		ll& dereference()const { return p_->data; }
		bool equal(const iterator& other) const { return p_ == other.p_; }
	};
	my_lsit(unique_ptr<my_list_node> p_) :p(move(p_)) {}
	iterator begin() { return iterator{ p.get() }; }
	iterator end() { return iterator(); }

private:
	unique_ptr<my_list_node> p;
};




#include <boost/iterator_adaptors.hpp>

class my_lsit2 {
public:
	// define iter 
	class iterator
		:public boost::iterator_adaptor<iterator, my_list_node*, ll, boost::forward_traversal_tag>
	{
	public:
		//iterator(my_list_node* f) : iterator_adaptor_(f) {}
		iterator() {}
		explicit iterator(my_list_node* p)
			: iterator::iterator_adaptor_(p) {}
		//ll& operator *() { return this->base()->data; }

	private:
		friend class boost::iterator_core_access;
		void increment() { this->base_reference() = this->base()->next.get(); }
		ll& dereference() const { return this->base()->data; }

	};
	my_lsit2(unique_ptr<my_list_node> p_) :p(move(p_)) {}
	iterator begin() { return iterator(p.get()); }
	iterator end() { return iterator(); }
private:
	unique_ptr<my_list_node> p;
};

my_lsit makelist(ll n) {
	unique_ptr<my_list_node> res;

	rep(i, 0, n) {
		res = make_unique<my_list_node>(i, move(res));
	}
	return my_lsit(move(res));
}
my_lsit2 makelist2(ll n) {
	unique_ptr<my_list_node> res;

	rep(i, 0, n) {
		res = make_unique<my_list_node>(i, move(res));
	}
	return my_lsit2(move(res));
}

int main_() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);

	auto list = makelist(10);

	for (auto x : list) {
		cout << x << endl;
	}

	auto list2 = makelist2(10);
	(list2.begin());
	for (auto x : list2) {
		cout << x << endl;
	}
	return 0;

}
