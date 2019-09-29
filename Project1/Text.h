#include "MyHeader.h"
#include "NumberTheory.h"

ll kmp_search(const string& text, const string& word) {
	// @return: the iterator where word was found in text
	if (word.empty()) return 0;

	vector<int> pi(word.size(), 0);
	for (int i = 1, k = 0; i < (int)word.size(); ++i) {
		while (k && word[k] != word[i]) k = pi[k - 1];
		if (word[k] == word[i]) ++k;
		pi[i] = k;
	}

	for (int i = 0, k = 0; i < (int)text.size(); ++i) {
		while (k && word[k] != text[i]) k = pi[k - 1];
		if (word[k] == text[i]) ++k;
		if (k == (int)word.size()) return i - k + 1;
	}
	return text.size();
}

template<ll MOD = MOD, ll B = 9973>
struct RollingHash {
	vector<modint<MOD>> hash, p;
	RollingHash() {}
	RollingHash(const string& s) :hash(s.size()+1), p(s.size()+1,1){
		rep (i ,0, s.size()) {
			hash[i + 1] = (hash[i] * B + s[i]);
			p[i + 1] = p[i] * B;
		}
	}
	//hash of S[l, r)
	modint<MOD> operator()(ll l, ll r) {
		return hash[r] - hash[l] * p[r - l];
	}
};
using RollingHash1 = RollingHash<MOD>;
using RollingHash2 = RollingHash<(1LL << 61) - 1>; 
using RollingHash3 = RollingHash<999999937>;

