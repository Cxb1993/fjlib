#ifndef FJLIBInterpolationH
#define FJLIBInterpolationH

#include "fjlib.h"
#include "fjlib_polynomial.h"
#include "fjlib_scale.h"

namespace fjlib {

/*!
// Parametric Interpolation Base Class
//
// Call set_data() to input the parametric scale(ordered) and data \n
// Call spline() to prepare if needed and then use splint() to 
// get interpolated value
//
//	1.23 revised
//	5.23 TFJScale_Mesh changed to TFJScale_GhostMesh
//	1.31.6 fixed a bug between(0.0-eps,1.0+eps) in splinx()
//
// Caution:
// splinx() is designed to capture roots for each segment
// but if the root happens on the node points, the result
// will be counted twice without further treatment
*/
class TFJInterp_Base {
public:
	typedef TFJScale_GhostMesh mesh_type;
	typedef TFJSeries series_type;
	typedef TFJPolynomial poly_type;
//	typedef std::pair<bool,bool> splinx_type;
protected:
	///	Internal vector of parametric scale
	vector_f		*_x;
	/// Uses TFJScale_Mesh to manage the parametric scale
	mesh_type		_scale;
	///	Internal vector of discretized data
	vector_f		*_y;
	/// Uses TFJSeries to manage the data value
	series_type		_series;
public:
	/// Set scale and data
	virtual
	void			set_data(vector_f* x,vector_f* y)
	{ 
		_x=x; _y=y; 
		_scale.link_vector(x); 
		_series.link_vector(y);
	}
	series_type		get_series() { return _series; }
	/// Get helper scale of the internal data
	mesh_type&		get_scale() { return _scale; }
	/// Call spline() to prepare the interpolation
	virtual
	void			spline()=0;
	/// Get interpolated value at scale=x
	virtual
	float_t			splint(float_t x)=0;
	/// Get interpolated value and derivative at scale=x
	virtual
	void			splint(float_t x, vector_f* v)=0;
	/// Returns the polynomial form for segment seg, using local coord [0,1]
	virtual
	poly_type		get_seg_poly(size_t seg)=0;
	/// Returns the scales at which values equal to y within segment seg
	/// Roots append to the vector instead of rewrite 
	virtual
	bool			splinx(size_t seg, float_t y, vector_f *v,
							float_t eps=0);
	/// Returns all the scales at which values equal to y
	virtual
	bool			splinx(float_t y, vector_f* v, float_t eps=0)
	{
/*
		v->resize(0);
		vector_f roots;
		for (size_t i=0; i<_series.seg_count(); i++)
		{
			splinx(i,y,&roots,eps);
			append(*v,roots);
		}
		return (v->size()>0);
*/
		v->resize(0);
		for (size_t i=0; i<_series.seg_count(); i++)
			splinx(i,y,v,eps);
		return (v->size()>0);
	}
};

/*!
// Interpolation class using polynomials for each segment
//
//	This is not a automated interpolation, you have to supply
//	both the data points and the polynomial which describe the 
//	segment and make sure the polynomial gives you same value
//	at the node point from two segment
//
//	Usage: the polynomial can be fed from other interpolation
//	class which ouputs polynomial, therefore it's possible to 
//	construct a mixed type interpolation scheme. 
//
*/
class TFJInterp_Poly: public TFJInterp_Base {
public:
	typedef std::vector<TFJPolynomial> poly_vec_type;
protected:
	/// Internal storage for all the polynomails
	poly_vec_type	_polys;
public:
	void			set_data(vector_f *x,vector_f* y)
	{ 
		TFJInterp_Base::set_data(x,y);
		_polys.resize(_scale.seg_count());
	}
    void			spline() {}
	float_t			splint(float_t x);
	void			splint(float_t x, vector_f* v);
	TFJPolynomial	get_seg_poly(size_t seg) 
	{ return _polys[seg]; }
	/// Access polynomial at each segment in local scale [0,1]
	virtual inline
	TFJPolynomial&	poly_at(size_t seg) { return _polys[seg]; }
	inline
	poly_vec_type&	polys() { return _polys; }
};


/*!
// Local Linear Interpolation Class
*/
class TFJInterp_Line: public TFJInterp_Base {
public:
    void			spline() {}		// no preparation here
	float_t			splint(float_t x);
	void			splint(float_t x, vector_f* v);
	TFJPolynomial	get_seg_poly(size_t seg);
//	bool			splinx(size_t seg, float_t y, vector_f *v);
//	bool			splinx(float_t y, vector_f* v) 
//	{ return TFJInterp_Base::splinx(y,v); }
};

}

//---------------------------------------------------------------------------
#endif
