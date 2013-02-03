#ifndef TFJLIBPolynomialH
#define TFJLIBPolynomialH

#include "fjlib.h"

namespace fjlib {

/*!
// Polynomial Class
//
// Support plus+, minus-, mutiply*, integrate<<, differential>>,
// evaluation|, shift%, scale& operations
// 
// Usage:
// You can call set_coef()  to setup coefficients or 
// call set_degree() and use operator[]() to access the coefficients.
// Operator +,-,*,>>,<<,|,&,% are supported.
//
// 1.21	revised, more operations added, template deleted
// 5.24 set_coef() for up to ten coefficients are added
*/
class TFJPolynomial {
public:
	typedef vector_f vec_type;
protected:
	/// Internal storage, vector based
	vec_type	vec;

	inline
	bool		outofbound(size_t deg)
	{ return (deg>=vec.size()); }

	/// [Helper] convert the size to degree
	inline
	size_t		deg2size(size_t d) const { return d+1; }
	/// [Helper] convert the degree to size
	inline
	size_t		size2deg(size_t s) const { return s-1; }
	void		set_zero_deg(float_t v)
	{
		vec.resize(1);	vec[0]=v;
	}
public:
	void		set_vector(const vector_f& v)
	{ vec=v; }

	/// Default constructor
	TFJPolynomial() { set_zero_deg(0.0); }
	/// Constructor, zero degree polynomial
	TFJPolynomial(float_t v) { set_zero_deg(v); }
	/// Constructor, assign
	TFJPolynomial(const TFJPolynomial& v)
	{ if (&v!=this) set_vector(v.vec); }
	/// Constructor for vector_f
	TFJPolynomial(const vector_f& v)
	{ if (&v!=&vec) set_vector(v); }
	/// Set coef
	void		set_coef(float_t v0)
	{ set_zero_deg(v0); }
	void		set_coef(float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1); }
	void		set_coef(float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2); }
	void		set_coef(float_t v3, float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2,v3); }
	void		set_coef(float_t v4, float_t v3, float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2,v3,v4); }
	void		set_coef(float_t v5, float_t v4, float_t v3, float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2,v3,v4,v5); }
	void		set_coef(float_t v6, float_t v5, float_t v4, float_t v3, float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2,v3,v4,v5,v6); }
	void		set_coef(float_t v7, float_t v6, float_t v5, float_t v4, float_t v3, float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2,v3,v4,v5,v6,v7); }
	void		set_coef(float_t v8, float_t v7, float_t v6, float_t v5, float_t v4, float_t v3, float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2,v3,v4,v5,v6,v7,v8); }
	void		set_coef(float_t v9, float_t v8, float_t v7, float_t v6, float_t v5, float_t v4, float_t v3, float_t v2, float_t v1, float_t v0)
	{ fjlib::set_vector(vec,v0,v1,v2,v3,v4,v5,v6,v7,v8,v9); }

public:
/* basic functions */
	/// Call set_degree() before using operator[]
	inline
	void		set_degree(size_t deg) 
	{ vec.resize(deg2size(deg)); }
	/// Returns the degree of the polynomial
	inline
	size_t		degree() const { return size2deg(vec.size()); }
	/// Returns the non-zero degree of the polynomial
	size_t		nz_degree() const;
	/// Returns the coefficient of degree deg, auto-expand supported
	inline
	float_t&	coef_at(size_t deg)
	{
		if (outofbound(deg)) vec.resize(deg2size(deg));
		return vec[deg];
	}
	/// Returns the coefficient of degree deg, auto_expand not supported but fast
	inline
	float_t&	operator[](size_t deg)
	{ return vec[deg]; }
	/// Returns the coefficients , const
	inline const
	float_t		operator[](size_t deg) const
	{ return vec[deg]; }
public:
/* hot functions */
	/// Evaluate the value at x
	float_t		eval(float_t x) const;
	/// eval(x1)-eval(x)
	float_t		eval_diff(float_t x0, float_t x1) const
	{ return eval(x1)-eval(x0); }
	/// Returns its 0-n derivatives value at x, operator | is supported
	void		eval(float_t x, size_t n, vector_f* v) const;
	/// Shift the abscissa by offset x<-x+offset
	void		shift(float_t offset);
	/// Scale the abscissa by factor x<-x*factor
	void		scale(float_t factor);
public:
/* operational functions */
	/// Operator for polynomial=polynomial
	TFJPolynomial&
				operator=(const TFJPolynomial& p)
	{ if (this!=&p) vec=p.vec; return *this; }
	/// Operator for polynomial+=polynomial
	TFJPolynomial&
				operator+=(const TFJPolynomial& rhs);
	/// Operator for polynomial+=value
	inline
	TFJPolynomial&
				operator+=(float_t rhs)
	{ vec[0]+=rhs; return *this; }
	/// Operator for polynomial*=polynomial
	TFJPolynomial&
				operator*=(const TFJPolynomial& rhs);
	/// Operator for polynomial*=value
	inline
	TFJPolynomial&        
				operator*=(float_t rhs)
	{ vec*=rhs; return *this; }
	/// Operator for polynomial>>=value, differentiation
	TFJPolynomial&        
				operator>>=(size_t rhs);
	/// Operator for polynomial<<=value, integration
	TFJPolynomial&
				operator<<=(size_t rhs);
	/// Operator for power(polynomial,k)
    TFJPolynomial&
				operator^=(size_t rhs);
	/// Operator for polynomial shift
	inline
    TFJPolynomial&
				operator%=(float_t rhs)
	{ shift(rhs); return *this; }
	/// Operator for polynomial scale
	inline
    TFJPolynomial&
				operator&=(float_t rhs)
	{ scale(rhs); return *this; }
};

std::ostream& operator<<(std::ostream& out, const TFJPolynomial& c);

TFJPolynomial operator+(const TFJPolynomial& rhs);
TFJPolynomial operator-(const TFJPolynomial& rhs);
TFJPolynomial operator+(const TFJPolynomial& lhs, const TFJPolynomial& rhs);
TFJPolynomial operator-(const TFJPolynomial& lhs, const TFJPolynomial& rhs);
TFJPolynomial operator*(const TFJPolynomial& lhs, const TFJPolynomial& rhs);
TFJPolynomial operator*(float_t lhs, const TFJPolynomial& rhs);
TFJPolynomial operator>>(const TFJPolynomial& lhs, size_t rhs);
TFJPolynomial operator<<(const TFJPolynomial& lhs, size_t rhs);
TFJPolynomial operator^(const TFJPolynomial& lhs, size_t rhs);
/// Operator for evaluating the polynomial at x=rhs
float_t operator|(const TFJPolynomial& lhs, float_t rhs);
TFJPolynomial operator%(const TFJPolynomial& lhs, float_t rhs);
TFJPolynomial operator&(const TFJPolynomial& lhs, float_t rhs);


typedef TFJPolynomial poly_f;
}

#endif

