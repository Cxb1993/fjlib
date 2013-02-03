#ifndef TFJLIBCubicSplineH
#define TFJLIBCubicSplineH

#include "fjlib_interpolation.h"

namespace fjlib {


/*! 
// Cubic Spline Class
//
// Call set_bc() to setup the boundary condition \n
// Cation: y" instead of y' is solved and stored
//
// 10.19 add set_lbc and set_rbc seperately
*/
class TFJCubicSpline: public TFJInterp_Base {
public:
	/// Boundary type for ends
	enum BCType {
		csbtNatural,	///< y''=Const, default is zero
		csbtClamped,	///< y'=Const
		csbtBessel,		///< y'=interpolated by three end points
		csbtQuadratic,	///< y'' same for two end points
		csbtClosed		///< y y' and y'' same for both ends
	};
private:
    vector_f		ap,al,ar,sr;
protected:
	vector_f		_y2;
	BCType			bt_b,bt_e;
	float_t			yp_b,yp_e;
public:
	TFJCubicSpline(): TFJInterp_Base(),
					bt_b(csbtNatural), bt_e(csbtNatural),
					yp_b(0), yp_e(0) {}
	void			set_bc(BCType b0,BCType b1,
							float_t v0, float_t v1)
	{ bt_b=b0; bt_e=b1; yp_b=v0; yp_e=v1; }
	void			spline();
	float_t			splint(float_t x);
	void			splint(float_t x, vector_f* v);
	TFJPolynomial	get_seg_poly(size_t seg);
/*
	bool			splinx(size_t seg, float_t y, vector_f *v,
							float_t eps=0);
	bool			splinx(float_t y, vector_f* v) 
	{ return TFJInterp_Base::splinx(y,v); }
*/

	///	Returns the internal y" value
	vector_f&		yp2() { return _y2; }
public:
	void			set_lbc(BCType b0, float_t v0)
	{ bt_b=b0; yp_b=v0; }
	void			set_rbc(BCType b1, float_t v1)
	{ bt_e=b1; yp_e=v1; }
};

}	// end of namespace

#endif

