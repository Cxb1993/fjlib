#ifndef FJLIB_INTERFACE_H
#define FJLIB_INTERFACE_H

#include "fjlib_curve.h"
#include "fjlib_scale.h"
#include "fjlib_smoothing.h"

namespace fjlib {

/*!
// Interface class
// A cubic spline class with remesh embeded
//
//	3.21.06 seperated from BubbleNS for reuse
//			removed axis related stuff
//			remove is_flat option
//			need to make apply_spline_bc() general later
*/
class TFJInterface {
public:
	typedef TFJCurve_CubicSpline	curve_type;
protected:
	curve_type		_curve;		// curve object
	vector_f		_ns;		// arc vector for new mesh
	vector_f		_xp1,_xp2;	// derivative of x(s)
	vector_f		_yp1,_yp2;	// derivative of y(s)
protected:	
	virtual void	initialize() { apply_spline_bc(); }
	virtual void	apply_spline_bc();
	virtual void	apply_bc();
	virtual	void	calc_derivs();
public:
	vector_f*		get_x() { return _curve.get_points().x(); }
	vector_f*		get_y() { return _curve.get_points().y(); }
	void			set_data(vector_f* x, vector_f* y)
	{ _curve.set_data(x,y); initialize(); }
	curve_type*		curve() { return &_curve; }
	/// spline the whole curve and its related properties
	virtual	void	spline();	
	/// get derivative based on curve arc length
	virtual	void	arc_deriv(const vector_f& v, size_t ord, 
							vector_f& ov);
	/// remesh is used to make a new mesh
	virtual	void	remesh_uniform(float_t npl);
	/// new mesh before finilized
	vector_f&		new_mesh_arc() { return _ns; }
	/// finalize_mesh copies new mesh over
	virtual	void	finalize_mesh() { _curve.remesh(_ns); }
public: // utilities
	/// returns node derivative vector x
	vector_f&		xp(size_t ord)
	{ 
		switch (ord) {
			case 1: return _xp1; break;
			case 2: return _xp2; break;
			default: throw; break;
		}
	}
	/// returns node derivative vector x
	vector_f&		yp(size_t ord)
	{
		switch (ord) {
			case 1: return _yp1; break;
			case 2: return _yp2; break;
			default: throw; break;
		}
	}
};

}	// end of namespace

#endif


