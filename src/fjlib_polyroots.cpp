#include "fjlib_polyroots.h"
#include <cmath>
#include <iostream>

namespace fjlib {

int TFJPolyRoots::solve(const TFJPolynomial& p)
{
	INF=false;
	rts.resize(0);
	size_t deg=p.nz_degree();
	switch (deg) {
		case 0: return solve_zero(p);
		case 1:	return solve_linear(p);
		case 2: return solve_quadratic(p);
		case 3: return solve_cubic(p);
		default: return solve_poly(p);
	};
}

vector_f& TFJPolyRoots::real_roots(vector_f& v)
{
	v.resize(0);
	if (rts.size()<1) return v;
	for (size_t i=0; i<rts.size(); i++)
		if (std::fabs(rts[i].imag())<eps) push_back(v,rts[i].real());
	return v;
}

int TFJPolyRoots::solve_zero(const TFJPolynomial &p)
{
	if (std::fabs(p[0])>eps) return 0;
	INF=true;
	return 0;
}

int TFJPolyRoots::solve_linear(const TFJPolynomial& p)
{
	add_root(-p[0]/p[1]);
	return 1;
}

int TFJPolyRoots::solve_quadratic(const TFJPolynomial& p)
{
	float_t a=p[2],b=p[1],c=p[0];
	cplx_f delta=cplx_f(sqr(b)-4.0*a*c);
	cplx_f q=(sqrt(delta)*sgn(b)+b)*cplx_f(-0.5);
	add_root(q/a);
	add_root(c/q);
	return 2;
}

int TFJPolyRoots::solve_cubic(const TFJPolynomial& p)
{
	float_t a=p[2]/p[3],b=p[1]/p[3],c=p[0]/p[3];
	float_t q=(sqr(a)-3.0*b)/9.0,
			r=(2.0*pow(a,3)-9.0*a*b+27.0*c)/54.0;
	float_t sg=sqr(r)-pow(q,3);
	if (sg<0)
	{
		float_t st=acos(r/pow(q,1.5));
		add_root(-2.0*sqrt(q)*cos(st/3.0)-a/3.0);
		add_root(-2.0*sqrt(q)*cos((st+2.0*M_PI)/3.0)-a/3.0);
		add_root(-2.0*sqrt(q)*cos((st-2.0*M_PI)/3.0)-a/3.0);
	}
	else
	{
		float_t aa=-sgn(r)*pow(fabs(r)+sqrt(sg),1.0/3.0);
		float_t bb=0;
		if (aa!=0) bb=q/aa;
		add_root(aa+bb-a/3.0);
		add_root(cplx_f(-0.5*(aa+bb)-a/3.0,sqrt(3.0)/2.0*(aa-bb)));
		add_root(cplx_f(-0.5*(aa+bb)-a/3.0,-sqrt(3.0)/2.0*(aa-bb)));
	}
	return 3;
}

}	// end of namespace
