#include<stdio.h>
#include<vector>
#include<algorithm>
using namespace std;
typedef long long ll;
typedef pair<ll, ll>pii;
typedef pair<ll, pii>pi3;
int dat[101010];
ll r = 1000000000000000000LL;
void isok(vector<pi3>v)
{
	fill(dat, dat + v.size(), 0);
	ll s = 0;
	for (int i = 0; i < v.size(); i++)
	{
		s += v[i].first;
		if (v[i].second.first == 0)dat[v[i].second.second] |= 1;
		else dat[v[i].second.second] |= 2;
	}
	int d[4];
	fill(d, d + 4, 0);
	for (int i = 0; i < v.size(); i++)d[dat[i]]++;
	if (d[3] != 0 || d[1] == 0 || d[2] == 0)r = min(r, s);
}
int main()
{
	struct A {
		struct AA{
			int AAa;
			string AAs;
		};
		struct BB {
			double BBb;
		};
	};
	A a;
	return 0;
}
