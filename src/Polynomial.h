#pragma once

#include "MyHeader.h"
#include "RecurrenceRelation.h"

template<class underlying_t = ll, ll max_dim = INF>
class SparsePolynomial
{
public:
	using coeff_t = vector<pair<ll, underlying_t>>;
	coeff_t coeff;

	static void contract(coeff_t& coeff) {
		if (coeff.empty()) return;
		sort_by<0>(coeff);
		auto cur = 0;
		rep(i,1,coeff.size()) {
			if (coeff[cur].first == coeff[i].first) {
				coeff[cur].second += coeff[i].second;
			}
			else {
				++cur;
				coeff[cur] = coeff[i];
			}
		}
		coeff.erase(coeff.begin() + cur + 1, coeff.end());
		reduce(coeff);
	}

	static void reduce(coeff_t& coeff) {
		if (coeff.empty()) return;
		auto cur = 0;
		rep(i, 0, coeff.size()) {
			if (coeff[i].second != 0) {
				coeff[cur] = coeff[i];
				++cur;
			}
		}
		coeff.erase(coeff.begin() + cur, coeff.end());
	}

	SparsePolynomial() : coeff() {}
	SparsePolynomial(int x) : coeff({ { 0, x } }) {}
	explicit SparsePolynomial(const coeff_t& coeff_) : coeff(coeff_) { sort_by<0>(coeff); }
	SparsePolynomial& operator<<=(ll d) {
		coeff_t res;
		for (auto it = coeff.rbegin(); it != coeff.rend(); ++it) {
			if (it->first + d <= max_dim) {
				res.push_back({ it->first + d, it->second });
			}
		}
		coeff = res;
		return *this;
	}
	friend SparsePolynomial operator<<(SparsePolynomial lhs, int rhs) {
		return lhs <<= rhs;
	}
	SparsePolynomial& operator*=(SparsePolynomial rhs) {
		// O(max(coeff.size() * (rhs.coeff.size())^2)
		coeff_t res;
		for(auto& i: coeff) {
			for(auto& j: rhs.coeff) {
				if (i.first + j.first <= max_dim)
					res.emplace_back(i.first + j.first, i.second * j.second);
			}
		}
		swap(coeff, res);
		contract(coeff);
		return *this;
	}
	friend SparsePolynomial operator*(SparsePolynomial lhs, const SparsePolynomial& rhs) {
		return lhs *= rhs;
	}
	SparsePolynomial& operator+=(const SparsePolynomial& rhs) {
		coeff.insert(coeff.end(), rhs.coeff.begin(), rhs.coeff.end());
		contract(coeff);
		return *this;
	}
	friend SparsePolynomial operator+(SparsePolynomial lhs, const SparsePolynomial& rhs) {
		return lhs += rhs;
	}
	bool operator==(const SparsePolynomial& rhs) const {
		return coeff == rhs.coeff;
	}
	underlying_t operator[](ll d) const {
		auto it = lower_bound(all(coeff), pair<ll, underlying_t>({ d, 0 }), [](auto&& l, auto&& r) { return l.first < r.first; });
		if (it != coeff.end() && it->first == d) {
			return it->second;
		}
		else {
			return 0;
		}
	}
};

template<class T = ll, ll max_dim = INF>
struct Polynomial {
	vector<T> a; // coefficients
	Polynomial() = default;
	Polynomial(const vector<T>& coeff) : a(coeff) {}
	Polynomial(ll c) : a{c} {}
	Polynomial(const vector<pair<ll, T>>& coeff) {
		//sort_by<0>(coeff);
		a.resize(coeff.back().first+1);
		rep(i, 0, coeff.size()) {
			a[coeff[i].first] = coeff[i].second;
		}
	}
	Polynomial& operator<<=(ll d) {
		auto n = a.size();
		a.resize(a.size() + d);
		rrep(i, 0, n) {
			a[i + d] = a[i];
		}
		rep(i, 0, d) {
			a[i] = 0;
		}
		return *this;
	}
	friend Polynomial operator<<(Polynomial lhs, int rhs) {
		return lhs <<= rhs;
	}
	Polynomial& operator>>=(ll d) {
		rep(i, 0, a.size() - d) {
			a[i] = a[i+d];
		}
		a.resize(a.size() - d);
		return *this;
	}
	friend Polynomial operator>>(Polynomial lhs, int rhs) {
		return lhs >>= rhs;
	}
	Polynomial& operator*=(Polynomial rhs) {
		auto res = vector<T>(min(SZ(a) + SZ(rhs.a)-1, max_dim+1));
		rep(i, 0, a.size()) {
			if (a[i] == 0)
				continue;
			rep(j, 0, rhs.a.size()) {
				if(i+j <= max_dim)
					res[i + j] += a[i] * rhs.a[j];
			}
		}
		a = res;
		return *this;
	}
	friend Polynomial operator*(Polynomial lhs, const Polynomial& rhs) {
		return lhs *= rhs;
	}
	Polynomial& operator+=(const Polynomial& rhs) {
		a.resize(max(a.size(), rhs.a.size()));
		rep(i, 0, rhs.a.size()) {
			a[i] += rhs.a[i];
		}
		return *this;
	}
	friend Polynomial operator+(Polynomial lhs, const Polynomial& rhs) {
		return lhs += rhs;
	}
	bool operator==(const Polynomial& rhs) const {
		return a == rhs.a;
	}
	function<T(ll)> div(const Polynomial& Q) const {
		// return generating function of rational function q(x)/p(x).
		// O((a.size())^2)
		
		auto q = a;
		auto p = Q.a;
		vector<T> tail;
		auto p0 = p[0];
		p.erase(p.begin());
		if (q.size() <= p.size()) {
			q.resize(p.size());
		}


		rep(i, 0, p.size()) {
			p[i] /= p0;
			p[i] *= -1;
		}
		rep(i, 0, q.size()) {
			q[i] /= p0;
		}
		auto& init = q;
		rep(i, 1, init.size()) {
			rep(j, 1, min(i,SZ(p))+1) {
				init[i] += p[j-1]* init[i-j];
			}
		}
		if (q.size() > p.size()) {
			tail.insert(tail.begin(), q.begin(), q.begin() + q.size() - p.size());
			q.erase(q.begin(), q.begin() + q.size() - p.size());
		}
		reverse(all(p));
		return [tail, init, p](ll n) {
			if (n < tail.size())
				return tail[n];
			else
				return solve_recurrence_relation(init, p, n - tail.size());
		};
	}
};


