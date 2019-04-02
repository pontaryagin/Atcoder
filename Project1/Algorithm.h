#pragma once
#include "MyHeader.h"

// Functor is expected to be functional<val ,bool>
// This function returns the maximum iterator that stisfies f(*it) == f(* inputs.begin())
template<
	typename Inputs,
	typename Functor,
	typename ValType = typename Inputs::value_type>
	pair<typename Inputs::iterator, typename Inputs::iterator>
	binary_solve(Inputs& inputs, Functor f)
{

	auto left = inputs.begin();
	auto right = inputs.end();
	auto n = inputs.size();

	auto left_val = f(*left);
	auto right_val = f(*(right - 1));

	// check 
	assert(n >= 2);
	assert(left_val != right_val);

	while (left + 1 != right)
	{
		auto mid = left + (right - left) / 2;
		if (f(*mid) == left_val)
		{
			left = mid;
		}
		else
		{
			right = mid;
		}
	}

	return { left,right };
}

template<typename T>
bool check_binary_solve(T left, T right) {
	return left != right;
}
template<>
bool check_binary_solve<double>(double left, double right) {
	return abs(left - right) > 0.0000000000001;
}

template<typename T, typename Functor>
pair<T, T > binary_solve(T left, T right, Functor f)
{
	auto left_val = f(left);
	auto right_val = f(right);

	// check 
	assert(left < right);
	if (right_val == left_val)
		throw invalid_argument("right_val == left_val");

	while (check_binary_solve<T>(left, right))
	{
		auto mid = left + (right - left) / 2;
		if (f(mid) == left_val)
		{
			left = mid;
		}
		else
		{
			right = mid;
		}
	}

	return { left,right };
}


template<int I>
vll proj(vpll v) {
	vll res(v.size());
	rep(i, 0, v.size()) {
		if (!I) res[i] = v[i].first;
		else res[i] = v[i].second;
	}
	return res;
}


template<int I, class T>
vll proj(T v) {
	vll res(v.size());
	rep(i, 0, v.size()) {
		res[i] = get<I>(v[i]);
	}
	return res;
}