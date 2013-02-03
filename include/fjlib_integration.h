#ifndef FJINTEGRATE_H
#define FJINTEGRATE_H

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "fjlib.h"
#include "fjlib_function.h"

namespace fjlib {

/*! 
//	Integration data structure
//
//	User defined data has to derive from this class
*/
struct TFJIntegrate_Data {
	/// Integration object
	void*		integs;	
	/// Set integration object
	void		set_obj(void *t) { integs=t; }
	/// Stores the coordinates for multi-dimension
	float_t		x[10];
};

/*!
//	1D Integration base class
//
//	Supports any user defined function(static or class function).
//	Default format is float(void*,float), but extra parameters or 
//	coefficients can be fed in by using the TFJIntegrate_Data and 
//	its derivative object.
*/
class TFJIntegrate_Base {
public:
	typedef TFJFunction1p func_type;
protected:
	func_type	function;
	TFJIntegrate_Data*
				data;
public:
	/// Dimension ID, reserved for multi-dimension integration
	int			dim_index;			// dimension index
	/// Set integral function
	virtual	inline
	void		set_function(func_type func)
	{ function=func; }
	/// Get the integral function
	virtual inline
	func_type	get_function() { return function; }
	/// Set integration data object
	virtual inline
	void		set_data(TFJIntegrate_Data* data_obj)
	{ data=data_obj; }
	/// Get integration data object
	virtual inline 
	TFJIntegrate_Data*
				get_data() { return data; }
	/// Default setup call
	virtual inline
	void		set(TFJIntegrate_Data* data_obj,
					func_type func, int dim_i=1)
	{
		set_function(func);
		set_data(data_obj);
		dim_index=dim_i;
	}
	/// Evaluate the function value at x
	virtual inline
	float_t		value_at(float_t* x)
	{ 
		data->x[dim_index-1]=*x;
		return function(data,x);
	}
	/// Integrate
	virtual float_t
				integrate()=0;
};

/*!
//	Integration class
//
//	Support integration from a to b, no rules defined yet
*/
class TFJIntegrate: public TFJIntegrate_Base {
public:
	/// Range of the integration
	float_t		a,b;
	/// Integrate from a  to b
	virtual float_t
				integrate(float_t x0, float_t x1)=0;
	float_t		integrate() { return integrate(a,b); }
};

/*! 
//	Trapzoid Integration class
//
*/
class TFJIntegrate_Trapzd: public TFJIntegrate {
public:
	/// Number of segments
	int			nseg;
	virtual float_t
				integrate(float_t x0, float_t x1, int n);
	float_t		integrate(float_t x0, float_t x1)
	{ return integrate(x0,x1,nseg); }
};

}	// end of namespace

#include "fjlib_gaussquad.h"

namespace fjlib {

/*!
//	Gauss Quadrature Integration class
*/
class TFJIntegrate_GaussQuad: public TFJIntegrate {
protected:
	///	Internal gauss weighting generator
	TFJGaussQuadrature
				gauss;		
public:
	///	Order of the GuassQuadrature
	int			order;
	TFJIntegrate_GaussQuad(): order(5) {}
	virtual float_t 
				integrate(float_t x0, float_t x1, int o);
	float_t		integrate(float_t x0, float_t x1)
	{ return integrate(x0,x1,order); }
};

}	// end of namespace

#endif