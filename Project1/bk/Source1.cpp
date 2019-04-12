#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
#include <cstdio>
#include <string>
#include <bitset>
#include <list>
#include <set>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <functional>
#include <queue>
// #include "bits/stdc++.h"

using namespace std;
//debug

#define rep(i, N) for (int i = 0; i < N; i++)
#define rep2(i, N, M) for (int i = N; i < M; i++)
#define pb push_back



typedef long long ll;
typedef pair<ll, ll> pll;
typedef pair<int, int> pii;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<string> vs;
typedef priority_queue<pll, vector<pll>, greater<pll>> greater_pq;


#define all(a)  (a).begin(),(a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define vec(a) vector<a>

#define BIT_NUM 32
ll f( ll x) {
	int cnt = 0;
	while (x  >= (170LL<<cnt)) {
		cnt++;
	}
	return 30 *(1LL<< cnt) + max((x >> cnt) - 140, 0LL) + 1;

}

int main_() {

	while (true) {
		ll c, d;
		cin >> c >> d;
		cout << f(d) - f(c)+1<<endl;
	}
	return 0;
}

