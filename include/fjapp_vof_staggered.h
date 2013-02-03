#ifndef TFJVOFGen2D_StaggeredH
#define TFJVOFGen2D_StaggeredH

#include "fjapp_vof.h"
#include "fjapp_staggered_axis.h"  

namespace fjlib {

/*!
// Staggered special case
// Use the mesh twice refined as the input one, and output
// the staggered vof values
// 2.21.06 add function fill_unknowns() accordingly from TFJVOFGen_Full
*/
template <class ST=TFJScale_Staggered_UniMesh,
			class CT=TFJCurve_CubicSpline>
class TFJVOFGen_Staggered {
public:
	typedef TFJVOFGen_Full<ST,CT>			vof_type;
	typedef typename vof_type::axis_type	axis_type;
	typedef typename vof_type::curve_type	curve_type;
	typedef typename vof_type::vof_curve_type	
											vof_curve_type;
	typedef TFJStaggeredPosType				staggered_type;
	typedef TFJStaggeredAxis_Helper<ST>		staggered_axis_helper;
protected:
	///	pointer to the axis object
	axis_type		*axis;
	///	Internal pointer to the curve object
	curve_type		*curve;
	/// Internal vof object
	vof_type		vof;
	/// Internal storage for the full vof at refined mesh
	matrix_f		mat;
	/// Internal refined axis 
	axis_type		naxis;
	/// Returns the cell vof and area based on the refined mesh
	virtual
	float_t			cell_area(size_t m, size_t n);
	/// Return the starting segments m,n for refined mesh, could be negtive
	void			staggered_mn(size_t m, size_t n, 
								staggered_type stag_loc,
								size_t cell_loc,
								int& nx, int& ny)
	{
		nx=(int)m*2-1, ny=(int)n*2-1;
		switch(stag_loc) {
			case sptN: break;
			case sptO: nx++; ny++; break;
			case sptU: ny--; break;
			case sptV: nx--; break;
			case sptP: nx--; ny--; break;
			default: throw; break;
		}
		switch(cell_loc) {
			case bptSW: break;
			case bptSE: nx++; break;
			case bptNW: ny++; break;
			case bptNE: nx++; ny++; break;
			default: throw; break;
		}
	}
public:
	/// Set the 2D axis and the curve
	void			set_data(axis_type* _axis, curve_type* _curve);
	void		 	vof_gen(bool fill=true);
 
	vof_type&		get_vof() { return vof; }
//	matrix_f&		get_vofmat() { return mat; }
	axis_type*		get_axis() { return axis; }
	curve_type*		get_curve() { return curve; }
	void			fill_unknowns(matrix_n &umat);

	// Get data using staggered notation, u,v,p,node
	float_t			get_staggered_vof(size_t m, size_t n, 
						staggered_type stag_loc);
	// Check if cell intersect with the curve
	bool			is_colored(size_t m, size_t n,
						staggered_type stag_loc,
						size_t cell_loc);
	/// Get bounding curve, call is_colored first, only one cell
	vof_curve_type&	vof_curve() { return vof.vof_curve(); }

	void			set_coordinates(bool use_cylinder_coordinates)
	{ vof.cylinder=use_cylinder_coordinates; }
private:
	void 			convert_xtoseg(const vector_f& v1, vector_n& v2); 
};

/*!  
//	Special case 
//	where surface force can be calcuated for each staggered cell
//	in cylindrical coordinates
//	6.25 variable surface tension added
//	6.26 segment detection algorithm improved to be general
//			, previous hack removed
*/
template <class ST=TFJScale_Staggered_UniMesh,
			class CT=TFJCurve_CubicSpline>
class TFJVOFGen_Staggered_ST:public TFJVOFGen_Staggered<ST,CT> {
protected:
	/// Calc UST for refined mesh
	virtual
	void			get_cell_st(size_t m,size_t n, 
								float_t& ust, float_t& vst);
	float_t			Bo;
protected:
	TFJInterp_Line*	_sigInterp;
	virtual
	void			get_cell_stg(size_t m,size_t n, 
								float_t& ust, float_t& vst);
public:
	void			set_Bo(float_t b) { Bo=b; }
	void			get_staggered_st(size_t m, size_t n,
		typename TFJVOFGen_Staggered<ST,CT>::staggered_type stag_loc,
						float_t& ust, float_t& vst);
	TFJVOFGen_Staggered_ST(): TFJVOFGen_Staggered<ST,CT>(),
							Bo(0), _sigInterp(NULL) {}
public:
	bool			uniform_st() { return (_sigInterp==NULL); }
	void			set_STInterp(TFJInterp_Line* interp)
	{ _sigInterp=interp; }
};

}	// end of namespace

#include "fjapp_vof_staggered.cpp" // for template issue

#endif

