#ifndef TFJVOFGen2DH
#define TFJVOFGen2DH

#include "fjlib_scale.h"
#include "fjlib_curve.h"
#include <map>

namespace fjlib {

/// Position class which defines the edges,corners and center of the box
enum TFJVOFBoxPosType {
	bptSW=0,bptS,bptSE,bptE,bptNE,bptN,bptNW,bptW,bptO
};

/*! 2D VOF Helper class
//	This application class is to generate curves to 
//	help calculate the volume of fraction of the cells
//	through Green's theory. \n
//	Basically, a TFJCurve intersects with a 2D grid TFJScales 
//	seperating the domain into internal and external regions, 
//	and then for each cell we can output a curve which
//	surrounds the internal(or external) part of the cell. \n
//	axis_type and curve_type are not templated but can be 
//	modified to any TFJCurve_Spline and TFJGroup<TFJScale> for now \n
//	eps value (1e-12, default) is used twice to compesate the trucation error:
//	1. root+(-)eps is used to guarantee intersets are sucessful within each cell
//	even at extreme situations, ex. when root are on the nodes of the scales
//	2. for the bounding curve, roots are all unique within eps to 
//	exclude the same root on nearby edges. 
//	The above turn out to be not contradict with each other  
//	but still cannot guaranteed that the bounding curve
//	only contains unique nodes due to that the final assembling
//	process doesn't exclude identical points.
//	5.24 template added, non-template version is at backup folder
//		now AxisType and CurveSpline type can vary
//		especially non-uniform mesh is supported now
//	5.25 changed template variable from axis to scale
//		it's more computer parsable
//	5.26 guarantee the uniqueness of roots in the same cell
//		fix a bug in the inside_box() caused by numerical trucation
//	5.27 fix a huge extreme bug, that only happens when two lines 
//		perpendicular with each other, where the tolerace calculation 
//		turns out to be wrong in _troot()
//		 fix a bug about ambiguous between calculated 0 and empty value
//		change empty value to be default -1.0 now
//	6.26 add arc_info as vof_curve() parameter, which records
//		original curve arc length of the bounding curve
*/
template <class ST=TFJScale_UniMesh,
			class CT=TFJCurve_CubicSpline>
class TFJVOFGen {
public:
	typedef ST							scale_type;
	typedef CT							curve_type;
	typedef TFJGroup_XYZ<ST>			axis_type;

//	typedef TFJAxis_UniMesh				axis_type;
//	typedef TFJCurve_Line				curve_type;	
//	typedef TFJCurve_CubicSpline		curve_type;	

	typedef TFJCurve_Poly				vof_curve_type;
	typedef std::map<size_t,vector_f>	seg_map_type;
	typedef TFJCurve::point_type		point_type;
	typedef TFJGroup_XYZ<vector_f>		pos_vec_type;
	
	typedef typename curve_type::interp_type
										interp_type;
protected:
	///	Internal pointer to the axis object
	axis_type		*axis;
	///	Internal pointer to the curve object
	curve_type		*curve;
	/// Internal mapping from segment id to intersects points
	seg_map_type	seg_map;
	/// Internal storage for ouput curve position
	pos_vec_type	ppos;
	/// Internal ouput curve object
	vof_curve_type	pcurve;

	/// Helper convert from box configure id to seg id
	size_t			stoid(size_t xseg, size_t yseg, size_t loc)
	{ 
		size_t n=xseg_count();
		size_t base=yseg*(2*n+1)+xseg;
		switch (loc) {
			case bptS: return base; break;
			case bptE: return base+n+1; break;
			case bptN: return base+2*n+1; break;
			case bptW: return base+n; break;
			default: throw; break;
		}
	}	
	/// Helper convert from box configure id to conner location
	point_type		ctopos(size_t xseg, size_t yseg, size_t loc)
	{
		point_type p(2);
		size_t sftx=0,sfty=0;
		switch (loc) {
			case bptSW: break;
			case bptSE: sftx++; break;
			case bptNE: sftx++; sfty++; break;
			case bptNW: sfty++; break;
			default: throw; break;
		};
        p.x()=axis->x().value_at(xseg+sftx);
		p.y()=axis->y().value_at(yseg+sfty);
		return p;
	}
private:
	size_t			ii,jj;
	vector_f		vrts; 
	vector_n		vloc;   
	// storage for original arc length in new vof curve
	vector_f		vof_arc;
private:
	/// root+-eps, take care of trucation error
	vector_n		_troot(float_t root, scale_type& scl);
	/// check if root is unique in a vector
	bool			_root_unique(float_t root, vector_f& vec);
	/// check if root is unique inside a cell
	bool			_root_unique(float_t root, int m, int n);
	/// push poly section value in, need hack to work
	void			_push_poly_value(size_t xy, const TFJPolynomial &p);
protected:	
//	vector_f		intercept(interp_type& itp, float_t v, float_t eps);
	///	Tags all the intersects points along the edge of the cell
	void			tag_box();
	///	Returns the postion of points along the cell
	point_type		box_pt(size_t i);
	///	Links all  the points together by inserting the segments
	void			link_box();
	/// Sort the points by their arc length
	void			sort_box();
	/// Inserts all conner points into intersects points
	void			corner_box();
	/// Helper to circle iterate through the box
	size_t			inc(size_t i);
	/// Insert a point & arc info into the output curve
	void			insert_pt(const point_type& pt,
							const float_t arc);
	/// Insert a line into the bounding curve
	void			insert_line(size_t in, size_t jn);
	/// Insert a polynomial into the bounding curve
	void			insert_poly(size_t in, size_t jn);
	/// Helper to check if any couple of intersects are inside the box
	bool			inside_box(size_t in, size_t jn);
	/// Helper to construct the mapping between edge and intersects
	void			insert_pair(size_t key, float_t v);
public:
	/// Set the 2D axis and the curve
	void			set_data(axis_type* _axis, curve_type* _curve)
	{
		axis=_axis; curve=_curve;
	}
	/// Generate edge and intersects mapping
	/*!
	//	Generate the mapping from seg id to roots, and
	//	seg id can be obtained using stoid()
	*/
	void			vof_gen();
	/// Returns if the box is colored
	bool			is_colored(size_t in, size_t jn);
	///	Outputs the bounding curve which is checked by is_colored()
	vof_curve_type& vof_curve(vector_f *arc_info=NULL);
	inline
	size_t			xseg_count() 
	{ return axis->at(0).seg_count(); }
	inline
	size_t			yseg_count() 
	{ return axis->at(1).seg_count(); }

	///	Ouputs the intersects at the segment (i,j,loc)
	vector_f		roots_at(size_t xseg, size_t yseg, size_t loc)
	{
		seg_map_type::iterator iter=seg_map.find(stoid(xseg,yseg,loc));
		if (iter!=seg_map.end())	// has roots on that seg
			return (*iter).second;
		else
			return vector_f();		
	}
public:
	/// Eps value which used for compensating the trucation error
	float_t			eps;
	/// value for cell which doesn't have a value
	float_t			void_value;
	TFJVOFGen(): eps(1e-12), void_value(-1.0) {}
	axis_type*		get_axis() { return axis; }
	curve_type*		get_curve() { return curve; }
protected:
	inline
	float_t			get_cell_area(int xseg, int yseg)
	{
		return axis->x().seg_len(xseg)*
			axis->y().seg_len(yseg);
	}
};
 

/*!
//	This class is not general and only specific for my problem
//	It can be made general by setting up callback function
//	Usage, input a matrix and output a full vof map 
*/
template <class CT=TFJCurve_CubicSpline,
			class AT=TFJAxis_UniMesh>
class TFJVOFGen_Full: public TFJVOFGen<CT,AT> {
protected:
	/// get unique sorted center roots at y=const
	void			get_ycenter_roots_xseg(size_t ys, vector_n& xrts);
public:
	virtual
	float_t			get_cell_vof(size_t m,size_t n);
	virtual
	void			calc_vofs(matrix_f &mat);
	virtual
	void			fill_vofs(matrix_f &mat);

	vector_f		test(size_t yseg)
	{
		vector_n a; get_ycenter_roots_xseg(yseg,a);
		return a; 
	}
	/// Use cylinder coordinates instead of cartisian
	bool			cylinder;

	TFJVOFGen_Full(): TFJVOFGen<CT,AT>(), cylinder(false) {}
	void			set_coordinates(bool use_cylinder_coordinates)
	{ cylinder=use_cylinder_coordinates; }
};
}	// end of namespace

#include "fjapp_vof.cpp"

#endif
