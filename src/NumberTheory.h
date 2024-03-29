#pragma once
#include "MyHeader.h"

// === Modint ===
//static uint_fast64_t runtime_modulus = MOD;

template <ll modulus>
class modint
{
public:
    ll val;
    constexpr modint() : val(0) {}
    constexpr modint(ll x) : val((x %= mod()) < 0 ? x + mod() : x) {}
    constexpr modint(ll x, ll modulus_) {
        set_modulo(modulus_); val = (x %= mod()) < 0 ? x + mod() : x;
    }
    template<class Ret = ll &>
    static auto modulo() -> std::enable_if_t<(modulus <= 0), Ret> {
        static ll runtime_modulus= numeric_limits<ll>::max(); return runtime_modulus; // singleton technique
    }
    template<class Ret = const ll>
    static auto mod() -> std::enable_if_t<(modulus <= 0), Ret> { return modulo(); }

    template<class Ret = const ll>
    static constexpr auto mod()->std::enable_if_t<(modulus > 0), Ret> { return modulus; }

    template<ll modulus_ = modulus, enable_if_t<(modulus_ <= 0), nullptr_t> = nullptr >
    static void set_modulo(ll mod) { modulo() = mod; }
    void reset_modulo(ll modulus_) { modulo() = modulus_; val %= mod();}
    constexpr modint inv() { return pow(mod() - 2); }
    constexpr ll value() const noexcept { return val; }
    constexpr modint operator+(const modint rhs) const noexcept {
        return modint(*this) += rhs;
    }
    constexpr modint operator-(const modint rhs) const noexcept {
        return modint(*this) -= rhs;
    }
    constexpr modint operator*(const modint rhs) const noexcept {
        return modint(*this) *= rhs;
    }
    constexpr modint operator/(const modint rhs) const noexcept {
        return modint(*this) /= rhs;
    }
    modint &operator+=(const modint rhs) noexcept {
        val += rhs.val;
        if (val >= mod()) {
            val -= mod();
        }
        return *this;
    }
    modint &operator-=(const modint rhs) noexcept {
        if (val < rhs.val) {
            val += mod();
        }
        val -= rhs.val;
        return *this;
    }
    modint &operator*=(const modint rhs) noexcept {
        val = val * rhs.val % mod();
        return *this;
    }
    modint &operator/=(modint rhs) noexcept {
        ll exp = mod() - 2;
        while (exp) {
            if (exp % 2) {
                *this *= rhs;
            }
            rhs *= rhs;
            exp /= 2;
        }
        return *this;
    }
    modint &operator++() noexcept {
        return *this += modint(1);
    }
    modint operator++(int) noexcept {
        modint t = *this;
        *this += modint(1);
        return t;
    }
    modint &operator--() noexcept {
        return *this -= modint(1);
    }
    modint operator--(int) noexcept {
        modint t = *this;
        *this -= modint(1);
        return t;
    }
    constexpr modint operator-() { return val ? mod() - val : val; }
    constexpr bool operator==(const modint rhs) const noexcept { return val == rhs.value(); }
    constexpr bool operator!=(const modint rhs)const  noexcept { return val != rhs.value(); }
    constexpr bool operator <(const modint rhs)const  noexcept { return val < rhs.value(); }
    static constexpr modint zero() { return modint(0); }
    static constexpr modint unit() { return modint(1); }

    modint pow(long long k) const {
        modint v = *this;
        modint res(1), tmp(v);
        while (k) {
            if (k & 1) res *= tmp;
            tmp *= tmp;
            k >>= 1;
        }
        return res;
    }

    ll log(modint b) {
        modint val = *this;
        const ll sq = 40000;
        map<modint, ll> dp;
        //dp.reserve(sq);
        modint res(1);
        for (ll r = 0; r < sq; r++) {
            if (!dp.count(res)) dp[res] = r;
            res *= val;
        }
        modint p = val.inv().pow(sq);
        res = b;
        for (ll q = 0; q <= mod() / sq + 1; q++) {
            if (dp.count(res)) {
                ll idx = q * sq + dp[res];
                if (idx > 0) return idx;
            }
            res *= p;
        }
        return INF;
    }
    friend ostream& operator <<(ostream& o, const modint<modulus>& t) {
        o << t.value();
        return o;
    }
    friend istream& operator >>(istream& in, modint<modulus>& t) {
        ll x;
        in >> x;
        t = modint<modulus>(x);
        return in;
    }
    friend modint<modulus> POW(modint<modulus> x, ll n) {
        return modint<modulus>(POW(x.value(), n, mod()));
    }


};

template<class T>
class Combination {
    // this calculates combination (nCk).
    // Constructor runs in O(MAX).
    // get(n,k) returns nCk in O(1).
    ll N_MAX;
    vector<T> fac, finv;
public:
    Combination(ll N_MAX = 210000)
        : N_MAX(2), fac(2,1), finv(2,1){
        pre_process(N_MAX);
    }
    T operator()(ll n, ll k) {
        // choose k from n
        pre_process(n);
        if (0<= n && n < k) return 0;
        if (k == 0) return 1;
        if (n < 0) return operator()(-n+k-1, k)* (k%2?-1:0);
        return fac[n] * (finv[k] * finv[n - k]);
    }
    T H(ll n, ll k) {
        // 1) 区間[0, k) を（空を許して）n個に分割する場合の数
        // 2) n個の中からk個を重複を許して選ぶ
        return operator()(n + k - 1, k);
    }
    T P(ll n, ll k) {// n (n-1) ... (n-k+1)
        pre_process(n);
        return (n<k|| n<0 )? T(0) : fac[n] * finv[n - k];
    }
    T Fac(ll n) { return P(n,n); }
    T FacInv(ll n) { pre_process(n);  return n < 0 ? T(0) : finv[n]; }
private:
    void pre_process(ll n) {
        if (n <= N_MAX) return;
        ll next = N_MAX;
        while (next <= n) next <<= 1;
        fac.resize(next); finv.resize(next);
        rep(i, N_MAX, next) {
            fac[i] = fac[i - 1] * i;
            finv[i] = finv[i - 1] / i;
        }
        N_MAX = next;
    }
};


ll choose(int n, int r) { // O(r) for small n
    ll acc = 1;
    rep(i, 0, r) acc = acc * (n - i) / (i + 1);
    return acc;
}

ll gcd(ll val, ll b) {
    if (val == INF) return b;
    if (b == INF) return val;
    if (val < 0) val *= -1;
    if (b < 0) b *= -1;
    if (val == 0) return b;
    if (b == 0) return val;
    if (val%b == 0) return b;
    else return gcd(b, val%b);
}




vll divisor(ll n) {
    // returns common divisors in O(sqrt(min(n,m)))
    vll res;
    for (ll i = 1; i*i <= n; i++) {
        if (i * i == n) res.push_back(i);
        else if (n%i == 0) {
            res.push_back(i);
            res.push_back(n / i);
        }
    }
    sort(res.begin(), res.end());
    return res;
}

vll divisor(ll n, ll m) {
    // returns common divisors in O(sqrt(min(n,m)))
    if (n > m) swap(n, m);
    vll res;
    for (ll i = 1; i*i <= n; i++) {
        if (n%i == 0) {
            if (m%i == 0) res.push_back(i);
            if (m % (n / i) == 0) res.push_back(n / i);
        }
        if (i * i == n && m % i == 0) res.push_back(i);
    }
    sort(res.begin(), res.end());
    return res;
}

vpll prime_factorize(ll n) {
    // returns prime factorization of n in O(sprt(n))
    vector<pll> res;
    for (ll p = 2; p*p <= n; ++p) {
        if (n%p != 0)continue;
        ll num = 0;
        while (n%p == 0) { ++num; n /= p; }
        res.push_back({ p,num });
    }
    if (n != 1) res.push_back(make_pair(n, 1));
    return res;
}

ll eulers_phi(ll n) {
    // returns the number of coprime integers in O(sqrt(n) + n^(1.4/log log n)) [ the number of prime divisors <= O(log(n^(1.4/log log n) ) ]
    ll res = 0;
    auto primes = prime_factorize(n);
    ll psize = primes.size();
    rep(bit, 1, 1LL << psize) {
        ll bitsize = popcnt(bit);
        ll arg = (bitsize % 2 ? 1 : -1);
        ll base = 1;
        rep(i, 0, primes.size()) {
            ll p = primes[i].first;
            if(bit &(1LL << i))
                base *= p;
        }
        res += arg * (n / base);
    }
    return n - res;
}

ll floor_div(ll a, ll b){
    int d = a / b;
    int r = a % b;
    return r ? (d - ((a < 0) ^ (b < 0))) : d;
}

ll ceil_div(ll a, ll b){
    int d = a / b;
    int r = a % b;
    return r ? (d + ((a < 0) ^ (b > 0))) : d;
}
