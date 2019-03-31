#pragma once
#include "MyHeader.h"


template<class T>
T atbit(T n, T i) {
	return (n >> i) % i;
}
template<class T>
T getbit(T i) {
	return 1LL << i;
}

