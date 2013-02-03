#ifndef TFJLIBPolyRootsH
#define TFJLIBPolyRootsH

// 5.24 solve_zero() added to cover extreme case

#include "fjlib.h"
#include "fjlib_polynomial.h"
#include <complex>

namespace fjlib {

class TFJPolyRoots {
public:
	typedef std::complex<float_t> cplx_f;
	typedef boost::numeric::ublas::vector<cplx_f> vec_cplx_f;
protected:
	vec_cplx_f		rts;
//	vector_f		real_rts;
	inline
	void			add_root(const cplx_f& v)
	{ push_back(rts,v); }
public:
	vec_cplx_f&		roots() { return rts; }
	virtual
	vector_f&		real_roots(vector_f& v);
	size_t			roots_count() { return rts.size(); }
	virtual
	int				solve(const TFJPolynomial& p);
protected:
	virtual
	int				solve_zero(const TFJPolynomial &p);
	virtual
	int				solve_linear(const TFJPolynomial &p);
	virtual
	int				solve_quadratic(const TFJPolynomial &p);
	virtual
	int				solve_cubic(const TFJPolynomial &p);
	virtual		
	int				solve_poly(const TFJPolynomial &p) { return 0; }
public:
	float_t			eps;
	TFJPolyRoots(): eps(1e-12), INF(false) {}
	bool			INF;
};

}	// end of namespace

#endif

