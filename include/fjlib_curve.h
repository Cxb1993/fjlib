#ifndef FJLIBCurveH
#define FJLIBCurveH

#include "fjlib_scale.h"
#include "fjlib_interpolation.h"

namespace fjlib {

/*! 
//	Curve Class
//	
//	6.2.2005 add calc_arcs() to spline()
*/
class TFJCurve {
public:
	typedef TFJGroup_XYZ<vector_f*> point_list_type;
	typedef TFJGroup_XYZ<float_t>	point_type;
	typedef TFJScale_GhostMesh		scale_type;
protected:
	/// Use a list of vector pointer to store x,y,z data
	/*!
	// This data access works best for me, 'cause i need
	// to constantly apply operations on vectors \n
	*/
	point_list_type	
				points;
	/// Internal storage for arc length
	vector_f	arcs_data;
	/// Helper storage for arc length
	scale_type
				arcs_scale;
	virtual
	void		setup() { calc_arcs(); }
public:
	TFJCurve() { arcs_scale.link_vector(&arcs_data); }

	/// Returns the dimenstion 
	const
	size_t		dim() const { return points.dim(); }
	/// Setup 1D data
	void		set_data(vector_f *x)
	{ points.redim(1); points.x()=x; setup(); }
	/// Setup 2D data
	void		set_data(vector_f *x, vector_f *y)
	{ points.redim(2); points.x()=x; points.y()=y; setup(); }
	/// Setup 3D data
	void		set_data(vector_f *x, vector_f *y, vector_f *z)
	{ points.redim(3); points.x()=x; points.y()=y; points.z()=z; setup(); }

	const 
	point_list_type&
				get_points() const { return points; }
	const
	point_type  pt_at(size_t i) const
	{
		point_type p(dim());
		for (size_t j=0; j<dim(); j++)
			p.at(j)=at(j,i);
		return p;
	}
	/// Access the coordinates of points at dim th dimension
	inline
	float_t&	at(size_t dim, size_t i)
	{ return (*points.at(dim))[i]; }
	inline const
	float_t&	at(size_t dim, size_t i) const
	{ return (*points.at(dim))[i]; }
	/// Access the x coordinates of points 
	inline
	float_t&	x(size_t i) { return (*points.x())[i]; }
	/// Access the y coordinates of points 
	inline
	float_t&	y(size_t i) { return (*points.y())[i]; }
	/// Access the z coordinates of points 
	inline
	float_t&	z(size_t i) { return (*points.z())[i]; }

	/// Returns the points number
	inline
	size_t		pt_count() { return points.x()->size(); }
	/// Returns the segment count
	inline
	size_t		seg_count() { return pt_count()-1; }
	/// Returns the arc length of the points, can be automatic or user defined
	inline
	scale_type&	arcs() { return arcs_scale; }
	/// Returns the total length of the curve
	inline
	float_t		len() { return arcs_scale.len(); }
	/// Returns the length of each segment
	inline
	float_t		len(size_t i) { return arcs_scale.seg_len(i); }
	/// Update the arc position of each node
	virtual
	void		calc_arcs();
};

std::ostream& operator<<(std::ostream& out, const TFJCurve& c);
std::istream& operator>>(std::istream& in, TFJCurve& c);

/*!
//	Spline Curve class
//	
//	Can plug in any spline algorithm from the template
//	Support the following for now,
//	
//	typedef TFJCurve_Spline<TFJInterp_Line> TFJCurve_Line;\n
//	typedef TFJCurve_Spline<TFJInterp_Poly> TFJCurve_Poly;\n
//	typedef TFJCurve_Spline<TFJCubicSpline> TFJCurve_CubicSpline;
*/
template <class Itp>
class TFJCurve_Spline: public TFJCurve {
public:
	typedef Itp					interp_type;
	typedef TFJGroup_XYZ<Itp>	interp_list_type;
protected:
	interp_list_type
				_interps;
	void		set_interp_data()
	{
		_interps.redim(dim());
		for (size_t i=0; i<dim(); i++)
			interp(i).set_data(&arcs_data,points.at(i));
	}
	void		setup() 
	{ 
		TFJCurve::setup(); set_interp_data();
	}
public:

	/// Returns the interpolation object for ith dim
	inline
	Itp&		interp(size_t i) { return _interps.at(i); }
	/// Prepare the spline
	virtual
	void		spline()
	{
		calc_arcs();
		for (size_t i=0; i<dim(); i++)
			interp(i).spline();
	}
	/// Interpolate the value s from the curve
	virtual
	point_type	splint(float_t s)
	{
		point_type p(dim());
		for (size_t i=0; i<dim(); i++)
			p.at(i)=interp(i).splint(s);
		return p;
	}
	/// Returns the interpolation list object
	interp_list_type&	
				interps() { return _interps; }
	/// Redistribute the nodes based on new arc length
	virtual
	void		remesh(const vector_f& arcs);
};

template <class Itp>
void TFJCurve_Spline<Itp>::remesh(const vector_f& arcs)
{
	size_t dm=dim(); 

	// create and fill in new pt location
	size_t nn=arcs.size();
	TFJGroup_XYZ<vector_f> mat(dm);
	for (size_t i=0; i<dm; i++)
		mat.at(i).resize(nn);
	for (size_t i=0; i<nn; i++)
	{
		point_type p=splint(arcs[i]);
		for (size_t j=0; j<dm; j++)
			mat.at(j)[i]=p.at(j);
	}

	// copy it back and setup
	for (size_t i=0; i<dm; i++)
		*points.at(i)=mat.at(i);
	setup();
}


typedef TFJCurve_Spline<TFJInterp_Line> TFJCurve_Line;
typedef TFJCurve_Spline<TFJInterp_Poly> TFJCurve_Poly;

}	// end of namespace

#include "fjlib_cubicspline.h"
namespace fjlib {
typedef TFJCurve_Spline<TFJCubicSpline> TFJCurve_CubicSpline;
}

#endif

