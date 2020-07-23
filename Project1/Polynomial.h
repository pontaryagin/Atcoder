#pragma once

#include "MyHeader.h"
#include "RecurrenceRelation.h"

template<class underlying_t = ll, ll dim = INF>
class SparsePolynomial
{
	using coeff_t = vector<pair<ll, underlying_t>>;
	coeff_t coeff;
public:
	SparsePolynomial() : coeff() {}
	SparsePolynomial(int x) : coeff(1, { 0, x }) {}
	SparsePolynomial(const coeff_t& coeff_) : coeff(coeff_) { sort(all(coeff)); }
	SparsePolynomial& operator<<=(ll d) {
		rrep(i, 0, coeff.size()) {
			if (d + coeff[i].first > dim) {
				coeff.pop_back();
			}
			else {
				coeff[i].first += d;
			}
		}
		return *this;
	}
	friend SparsePolynomial operator<<(SparsePolynomial lhs, int rhs) {
		return lhs <<= rhs;
	}
	SparsePolynomial& operator*=(SparsePolynomial rhs) {
		// O(max(coeff.size() * (rhs.coeff.size())^2)
		SparsePolynomial res;
		rep(i, 0, rhs.coeff.size()) {
			ll cur = 0;
			ll res_org_size = res.coeff.size();
			rep(j, 0, coeff.size()) {
				ll out_d = rhs.coeff[i].first + coeff[j].first;
				if (out_d <= dim) {
					while (cur < res_org_size && res.coeff[cur].first < out_d) {
						++cur;
					}
					if(cur < res_org_size && res.coeff[cur].first == out_d){
						res.coeff[cur].second += rhs.coeff[i].second * coeff[j].second;
					}
					else {
						res.coeff.push_back({ out_d, rhs.coeff[i].second * coeff[j].second });
					}
				}
			}
			sort(all(res.coeff));
		}
		swap(res, *this);
		return *this;
	}
	friend SparsePolynomial operator*(SparsePolynomial lhs, const SparsePolynomial& rhs) {
		return lhs *= rhs;
	}
	SparsePolynomial& operator+=(const SparsePolynomial& rhs) {
		ll cur = 0;
		ll org_size = coeff.size();
		bool need_sort = false;
		rep(i, 0, rhs.coeff.size()) {
			while (cur < org_size && coeff[cur].first < rhs.coeff[i].first) {
				++cur;
			}
			if (cur >= org_size) {
				coeff.insert(coeff.end(), rhs.coeff.begin() + i, rhs.coeff.end());
				break;
			}
			else {
				if (coeff[cur].first == rhs.coeff[i].first) {
					coeff[cur].second += rhs.coeff[i].second;
				}
				else {
					coeff.push_back(rhs.coeff[i]);
				}
			}
		}
		if (need_sort) {
			sort(all(coeff));
		}
		return *this;
	}
	friend SparsePolynomial operator+(SparsePolynomial lhs, const SparsePolynomial& rhs) {
		return lhs += rhs;
	}
	bool operator==(const SparsePolynomial& rhs) const {
		return coeff == rhs.coeff;
	}
	underlying_t operator[](ll d) const {
		auto it = lower_bound(all(coeff), pair<ll, underlying_t>{ d, 0 },
			[](const pair<ll, underlying_t>& l, const pair<ll, underlying_t>& r) {return l.first < r.first; });
		if (it != coeff.end() && it->first == d) {
			return it->second;
		}
		else {
			return 0;
		}
	}
};

template<class T = ll>
struct Polynomial {
	vector<T> a; // coefficients
	Polynomial() = default;
	Polynomial(const vector<T>& coeff): a(coeff) {}
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
		auto res = vector<T>(a.size() + rhs.a.size()-1);
		rep(i, 0, a.size()) {
			rep(j, 0, rhs.a.size()) {
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


