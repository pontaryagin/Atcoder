// ConsoleApplication1.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

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
//#include <bits/stdc++.h>

using namespace std;
//debug
#define dump(x)  cerr << #x << " = " << (x) << endl;
#define debug(x) cerr << #x << " = " << (x) << " (L" << __LINE__ << ")" << " " << __FILE__ << endl;
int main() {

	int N;
	string S;
	cin >> N >> S;
	char c = S[0];
	int val =N-1- count(S.begin(), S.end(), c);
	int ret = val;
	for (int i = 1; i < N; i++) {
		if (S[i - 1] = S[i]) {
			if (S[i] == 'W') {
				val += 1;
			}
			else {
				val -= 1;
			}
		}
		ret = min(ret, val);
	}
	cout << ret;
	return 0; 
}