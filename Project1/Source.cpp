
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
#include <assert.h>

using namespace std;
//debug

#define rep(i, N) for (int i = 0; i < N; i++)

#define pb push_back

typedef long long ll;
typedef pair<ll, ll> pll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<string> vs;


#define all(a)  (a).begin(),(a).end()
#define rall(a) (a).rbegin(), (a).rend()



int main() {
	ll n,m,d;
	cin >> n >> m >> d;
	printf("%.8f", double(m - 1) / n / n*(n*ll(d == 0) + 2*(n - d)*ll(d != 0)));
	cin>>n;

}

