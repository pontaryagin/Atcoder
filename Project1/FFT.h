#include "NumberTheory.h"
#include <complex>

vll plus_v(const vll& rhs, const vll& lhs) {
	assert(rhs.size() == lhs.size());
	vll res(rhs.size());
	rep(i, 0, rhs.size()) {
		res[i] = rhs[i] + lhs[i];
	}
	return res;
}


namespace FFT {
	// Fast Fourier Transformation
	using F= double;
	using C= complex<F>;
	// calc descrete Fourier transform
	vll trasform(const vll& f) {
		
	}
	
	vll inverse(const vll& F) {

		
	}
	
	
}


namespace NTT {
	// Number Theoretical Transformation
	

}