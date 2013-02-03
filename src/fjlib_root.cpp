#include "fjlib_root.h"

namespace fjlib {

bool TFJRoot_Bracket::brac(float_t& _x1,float_t& _x2)
{
	const int NTRY=50;
	const float_t FACTOR=1.6;
	float_t f1,f2;

	if (_x1 == _x2) throw "\nBad initial range in fjlib_root.cpp\n";
	f1=call_xf(_x1);
	f2=call_xf(_x2);
	for (size_t j=0;j<NTRY;j++) {
		if (f1*f2 < 0.0) return true;
		if (std::fabs(f1) < std::fabs(f2))
			f1=call_xf(_x1 += FACTOR*(_x1-_x2));
		else
			f2=call_xf(_x2 += FACTOR*(_x2-_x1));
	}
	return false;
}

float_t TFJRoot_Bracket::bisect(float_t _x1,float_t _x2)
{
	const int JMAX=40;
	float_t dx,f,fmid,xmid,rtb;

	f=call_xf(_x1);
	fmid=call_xf(_x2);
	if (f*fmid >= 0.0) throw "\nRoot must be bracketed for bisection in fjlib_root.cpp\n";
	rtb = f < 0.0 ? (dx=_x2-_x1,_x1) : (dx=_x1-_x2,_x2);
	for (size_t j=0;j<JMAX;j++) {
		fmid=call_xf(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (std::fabs(dx) < m_fEps || fmid == 0.0) return rtb;
	}
	throw "\nToo many bisections in fjlib_root.cpp\n";
}

float_t TFJRoot_Brent::brent()
{
    using namespace std;
	const int ITMAX=100;
	const float_t EPS=std::numeric_limits<float_t>::epsilon();
	float_t a=m_fMin,b=m_fMax,c=m_fMax,d,e,min1,min2;
	float_t fa=call_xf(a),fb=call_xf(b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		throw "\nRoot must be bracketed in fjlib_root.cpp";
	fc=fb;
	for (size_t iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*m_fEps;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=call_xf(b);
	}
	throw "\nMaximum number of iterations exceeded in fjlib_root.cpp";
}

float_t TFJRoot_Newton::newton(float_t x1, float_t x2)
{
	const int JMAX=20;
	float_t df,dx,f,rtn;

	vector_f val;
	rtn=0.5*(x1+x2);
	for (size_t j=0;j<JMAX;j++) {
		f=call(rtn,&val);
		df=val[1];
//		_funcd(_obj,rtn,f,df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			throw "\nJumped out of brackets in fjlib_root.cpp\n";
		if (std::fabs(dx) < m_fEps) return rtn;
	}
	throw "\nMaximum number of iterations exceeded in fjlib_root.cpp\n";
}

float_t	TFJRoot_Newton::safe_newton(float_t x1, float_t x2)
{
    using namespace std;
	const int MAXIT=100;
	float_t df,dx,dxold,f,fh,fl,temp,xh,xl,rts;

//	float_t x1=low_bound,x2=upp_bound;
	vector_f val;
//	_funcd(_obj,x1,fl,df);
//	_funcd(_obj,x2,fh,df);
	fl=call_xf(x1);
	fh=call_xf(x2);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		throw std::out_of_range("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
//	_funcd(_obj,rts,f,df);
	call_xfdf(rts,f,df);
	for (size_t j=0;j<MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < m_fEps) return rts;
//		_funcd(_obj,rts,f,df);
		call_xfdf(rts,f,df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	throw std::out_of_range("Maximum number of iterations exceeded in rtsafe");
}

}	// end of namespace