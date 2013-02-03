#include "fjlib_polynomial.h"

namespace fjlib {

size_t TFJPolynomial::nz_degree() const
{
	int i=degree();
	while (i>=0)
	{
		if (vec[i]!=0) return i;
		i--;
	}
	return 0;
}

float_t	TFJPolynomial::eval(float_t x) const
{ // copied from NR c++ book
	size_t j=nz_degree();
	float_t eval=vec[j];
	while (j>0) eval=eval*x+vec[--j];
	return eval;
}

void TFJPolynomial::eval(float_t x, size_t n, vector_f* v) const
{  // copied from Numerical Recipes, ddpoly()
	int nnd,j,i;
	float_t cnst=1;

	if (v->size()<n) v->resize(n);
	int nc=nz_degree();
	int nd=n-1;
	(*v)[0]=vec[nc];
	for (j=1;j<nd+1;j++) (*v)[j]=0;
	for (i=nc-1;i>=0;i--) {
		nnd=(nd < (nc-i) ? nd : nc-i);
		for (j=nnd;j>0;j--)
			(*v)[j]=(*v)[j]*x+(*v)[j-1];
		(*v)[0]=(*v)[0]*x+vec[i];
	}
	for (i=2;i<nd+1;i++) {
		cnst *= i;
		(*v)[i]*= cnst;
	}
}

void TFJPolynomial::shift(float_t offset)
{ // copy from dealII package
	size_t nz=nz_degree();
	if (nz<1) return;
    vec_type nc=vec;
    for (size_t d=1; d<deg2size(nz); ++d)
	{
		const size_t n = d;
        size_t bc = 1;
        float_t offset_power = offset;
        for (size_t k=0;k<d;++k)
		{
			bc=(bc*(n-k))/(k+1);
            nc[d-k-1]+=nc[d]*bc*offset_power;
            offset_power*=offset;
        }
//		if (bc==1)
//			throw std::length_error("bc=1");
	}
	vec=nc;
}

void TFJPolynomial::scale(float_t factor)
{
    float_t f = 1;
	size_t i=0;
	while (i<deg2size(nz_degree()))
	{
		vec[i]*=f;
        f*=factor;
		i++;
	}  
}

TFJPolynomial& TFJPolynomial::operator+=(const TFJPolynomial& rhs)
{
	size_t osz=nz_degree(),
			rsz=rhs.nz_degree();
	if (osz<=rsz)
	{
		set_degree(rsz); 
		vec+=rhs.vec;
	}
	else
		for (size_t i=0;i<=rsz; i++)
			vec[i]+=rhs.vec[i];
	return *this;
}

TFJPolynomial& TFJPolynomial::operator*=(const TFJPolynomial& rhs)
{
	size_t osz=nz_degree(),
			rsz=rhs.nz_degree();
	vec_type vm(deg2size(osz+rsz));
	vm*=0.0;	// initial
	
	for (size_t i=0; i<=osz; i++)
		for (size_t j=0; j<=rsz; j++)
			vm[i+j]+=vec[i]*rhs.vec[j];
	vec=vm;
	return *this;
}

TFJPolynomial& TFJPolynomial::operator>>=(size_t rhs)
{
	if (rhs==0) return *this;
	size_t nz=nz_degree();

	if (nz<rhs) { set_zero_deg(0.0); return *this; } 

	for (size_t i=0; i<rhs; i++)
	{
		for (size_t j=0; j<nz; j++)
			vec[j]=vec[j+1]*(j+1);
		nz--;
	}
	set_degree(nz);
	return *this;
}

TFJPolynomial& TFJPolynomial::operator<<=(size_t rhs)
{
	if (rhs==0) return *this;
	size_t nz=nz_degree();

	set_degree(nz+rhs);
	for (size_t i=0; i<rhs; i++)
	{
		for (size_t j=nz+1; j>0; j--)
			vec[j]=vec[j-1]/j;
		vec[0]=0.0;
		nz++;
	}
	return *this;
}

TFJPolynomial& TFJPolynomial::operator^=(size_t rhs)
{
	if (rhs==0) { set_zero_deg(1.0); return *this; }
	if (rhs==1) return *this;
	for (size_t i=0; i<rhs-1; i++)
		(*this)*=(*this);
	return *this;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
TFJPolynomial operator+(const TFJPolynomial& rhs)
{
	return rhs;
}

TFJPolynomial operator-(const TFJPolynomial& rhs)
{
	TFJPolynomial b=rhs;
	for (size_t i=0; i<=b.nz_degree(); i++)
		b[i]=-b[i];
	return b;
}

TFJPolynomial operator+(const TFJPolynomial& lhs, const TFJPolynomial& rhs)
{
	TFJPolynomial b=lhs;
	b+=rhs;
	return b;
}

TFJPolynomial operator-(const TFJPolynomial& lhs, const TFJPolynomial& rhs)
{
	return lhs+(-rhs);
}

TFJPolynomial operator*(const TFJPolynomial& lhs, const TFJPolynomial& rhs)
{
	TFJPolynomial b=lhs;
	b*=rhs;
	return b;
}

TFJPolynomial operator*(float_t lhs, const TFJPolynomial& rhs)
{
	return rhs*lhs;
}

TFJPolynomial operator>>(const TFJPolynomial& lhs, size_t rhs)
{
	TFJPolynomial b=lhs;
	b>>=rhs;
	return b;
}

TFJPolynomial operator<<(const TFJPolynomial& lhs, size_t rhs)
{
	TFJPolynomial b=lhs;
	b<<=rhs;
	return b;
}

TFJPolynomial operator^(const TFJPolynomial& lhs, size_t rhs)
{
	TFJPolynomial b=lhs;
	b^=rhs;
	return b;
}

TFJPolynomial operator%(const TFJPolynomial& lhs, float_t rhs)
{
	TFJPolynomial b=lhs;
	b%=rhs;
	return b;
}

TFJPolynomial operator&(const TFJPolynomial& lhs, float_t rhs)
{
	TFJPolynomial b=lhs;
	b&=rhs;
	return b;
}

float_t operator|(const TFJPolynomial& lhs, float_t rhs)
{
	return lhs.eval(rhs);
}


std::ostream& operator<<(std::ostream& out, const TFJPolynomial& c)
{
	int i=c.degree();
	while (i>0)
	{
		if (c[i]>0) out<< '+';
		if (c[i]!=0)
			out << c[i] << 'x' << '^' << i;
		i--;
	}
	if (c[0]>0) out << '+' ;
	if (c[0]!=0) out << c[0]; 
	return out;
}


}
