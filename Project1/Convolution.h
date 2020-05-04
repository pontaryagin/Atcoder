#include "MyHeader.h"


namespace convolution {
	namespace gcd {
		// (a[i]) => (f(a)[i] = \sum_{k|i} a[i])
		template<class T>
		void trans(vector<T>& a) {
			rep(k, 1, a.size()) {
				for (ll i = k * 2; i < a.size(); i += k) {
					a[k] += a[i];
				}
			}
		}
		template<class T>
		void inv(vector<T>& a) {
			rrep(k, 1, a.size()) {
				for (ll i = k * 2; i < a.size(); i += k) {
					a[k] -= a[i];
				}
			}
		}
		template<class T>
		vector<T> convolute(vector<T> a, vector<T> b) {
			trans(a); trans(b);
			rep(i, 0, a.size()) {
				a[i] *= b[i];
			}
			inv(a);
			return a;
		}
	}
}
 