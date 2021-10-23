#pragma once

#include "MyHeader.h"

template<class T>
T solve_recurrence_relation(const vector<T>& a, const vector<T>& p, ll n) {
    // solve : a[n] = p[0] * a[n-k] + ... + p[k-1] * a[n-1] by calc coefficient (Kitamasa method)
    // a[n] = x[0] * a[0] + x[1] * a[1] + ... + x[k-1] * a[k]

    assert(a.size() == p.size());
    auto k = SZ(p);
    if (n < k) {
        return a[n];
    }
    auto increment = [&](vector<T>& x) {
        x.insert(x.begin(), 0);
        rep(i, 0, k) {
            x[i] += x.back() * p[i];
        }
        x.pop_back();
    };
    auto dbl = [&](vector<T>& x) {
        vector<T> res(k);
        auto x_i = x;
        vector<vector<T>> X(k);
        X[0] = x_i;
        rep(i, 1, k) {
            increment(x_i);
            X[i] = x_i;
        }
        rep(i, 0, k) {
            rep(j, 0, k) {
                res[i] += x[j] * X[j][i];
            }
        }
        swap(res, x);
    };
    function<void(vector<T>&, ll)> rec = [&](vector<T>& x, ll n) {
        if (n == 0) {
            x.resize(k);
            x[0] = 1;
        }
        else if (n & 1) {
            rec(x, --n);
            increment(x);
        }
        else {
            n /= 2;
            rec(x, n);
            dbl(x);
        }
    };
    vector<T> x;
    rec(x, n);
    T res = 0;
    rep(i, 0, k) {
        res += x[i] * a[i];
    }
    return res;
};
