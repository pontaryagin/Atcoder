#pragma once
#include "MyHeader.h"


template<ll bit = 2LL>
ll at_bit(ll n, ll i) {
	return n / POW(bit, i) % bit;
}

template<>
ll at_bit<2>(ll n, ll i) {
	return (n >> i) % 2LL;
}

template<ll bit>
ll get_bit(ll i) {
	return POW(bit, i);
}

template<>
ll get_bit<2>(ll i) {
	return 1LL << i;
}

template<ll bit = 2>
ll get_max_bit(ll n) {
	ll tmp = bit;
	ll at = 0;
	while (tmp <= n) {
		at++;
		tmp *= bit;
	}
	return at;
}

template<>
ll get_max_bit<2>(ll n) {
	ll tmp = 2;
	ll at = 0;
	while (tmp <= n) {
		at++;
		tmp <<= 1;
	}
	return at;
}