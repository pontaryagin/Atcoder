#include "MyHeader.h"

template<ll dim = INF>
class SparsePolynomial
{
	vpll coeff;
public:
	SparsePolynomial() : coeff() {}
	SparsePolynomial(int x) : coeff(1, { 0, x }) {}
	SparsePolynomial(const vpll& coeff_) : coeff(coeff_) { sort(all(coeff)); }
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
	SparsePolynomial& operator*=(const SparsePolynomial& rhs) {
		SparsePolynomial res;
		rep(i, 0, coeff.size()) {
			ll cur = 0;
			ll res_org_size = res.coeff.size();
			rep(j, 0, rhs.coeff.size()) {
				ll out_d = coeff[i].first + rhs.coeff[j].first;
				if (out_d <= dim) {
					while (cur < res_org_size && res.coeff[cur].first < out_d) {
						++cur;
					}
					if(cur < res_org_size && res.coeff[cur].first == out_d){
						res.coeff[cur].second += coeff[i].second * rhs.coeff[j].second;
					}
					else {
						res.coeff.push_back({ out_d, coeff[i].second * rhs.coeff[j].second });
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
				coeff.insert(coeff.end(), rhs.coeff.begin()+i, rhs.coeff.end());
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
};

