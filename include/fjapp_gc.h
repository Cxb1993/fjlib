#ifndef FJAPP_GC_IMMERSED
#define FJAPP_GC_IMMERSED

#include "fjlib_vecmat.h"
#include "fjapp_vof.h"

namespace fjlib {

/*!
// need to find a way to record nodes around the ghost node
!*/
struct TFJGCNodeInfo {
	float_t x,y;	// position
	float dist;		// distance from ghost node
	int	type;		// 0:any; 1:on curve; 2: on grid
	int i,j;		// index depends on type
					// n/a(0); curve index(1); grid index(2)
	float_t v;		// scale value
};

/*!
//	2D Ghost-cell immersed boundary treatment for complex geometry 
//	using a rectangular mesh and a spline curve,  base class
//	
// 	created 1/4/2006
//
//	Notes:
//		dimension actually is embedded inside template classes,
//	therefore scale is used instead axis type.
//	1.26.06 fixed a bug where zero area of interpolation
//	1.27.06 find a huge bug due to confusion about npos_type & pos_type
// 	1.30.06 resolved the above bug, to ease the confusion,
//			position of ghost and surface are added
//	2.6.06	add symmetry support on x
//	2.27.06 fix a bug related to empty ghost cell
//			change the order of ghost search by y instead of x
//	2.28.06 fix a bug of finding mirror node, reset eps be 0.1
!*/
template <class ST=TFJScale_UniMesh,
			class CT=TFJCurve_CubicSpline>
class TFJGCGen {
public:
	typedef ST							scale_type;
	typedef CT							curve_type;
	typedef TFJGroup_XYZ<ST>			axis_type;
	typedef TFJGCNodeInfo				node_type;
	typedef std::vector<TFJGCNodeInfo>	node_vec_type;
	typedef std::vector<node_vec_type>	node_mat_type;
	typedef TFJGroup_XYZ<float_t>		pos_type;
	typedef std::vector<pos_type>		pos_vec_type;
	typedef TFJGroup_XYZ<int>			npos_type;
	typedef std::vector<npos_type>		npos_vec_type;			
private:
	size_t			_bound_min(size_t v, size_t dv);
	/// cation, here m is size+1
	size_t			_bound_max(size_t v, size_t m);
	float_t 		_dist(float_t x1, float_t y1,
						float_t x2, float_t y2);
	bool 			_get_range(float_t x, float_t y, int n,
								size_t& mx1, size_t& mx2, 
								size_t& my1, size_t& my2);
	bool 			_sort_vec(const vector_f& v, size_t n,
								vector_n& od);
//	void			_get_ghosts_pos(pos_vec_type& vp);
protected:
	/// Internal axis object
	axis_type 		*_axis;
	curve_type		*_curve;
	matrix_n		_umat;
	/// Stores the coordinates index of the boundary nodes
	npos_vec_type	_bnodes; 
	/// Stores the coordinates index of the ghosts point
	npos_vec_type	_ghosts;
	/// Outline the unknowns map with bnodes and ghosts
	/// Also create _bnodes and _ghosts for storage
	void			gen_unknown_map();
	/// Ghost nodes envirment includes all the info of nearby nodes 
	node_mat_type	*_env_info;
	/// mirror surface node pos storage
	node_vec_type	_mirrors;		
	/// Create nodes that perbenticular to the surface and 
	void			gen_mirrors();
//	/// Create mirror ghost nodes which reside inside known domain
//	void			gen_mirror_ghosts();
public: 
	/// Search n nearest nodes from (x,y) on bulk and sppend to nvec
	void			search_grid_nodes(float_t x, float_t y, 
								node_vec_type& nvec, size_t n);
	/// Search n nearest nodes from (x,y) on curve and append to nvec
	void			search_surf_nodes(float_t x, float_t y,
								node_vec_type& nvec, size_t n);
public:	// for testing
	/// Set enviorment info storage pointer
	void			set_env(node_mat_type *env)
	{ _env_info=env; }
	/// Detect and outline all the nearby nodes of each ghost node
	void			gen_env_info(const pos_vec_type& vp,
									size_t bn, size_t sn,bool sort=false);
	/// Based on type, detect ghost or surf nodes
//	void			gen_env_info(size_t type, size_t bn, size_t sn);
	/// Fill scale values to eviorment nodes
	void			fill_env_nodes(const matrix_f& vb,
									const vector_f& vs);
	/// Gen interpolation coefficient
	void			gen_interp_coef(const vector_n& src,
									matrix_f& coef);
	/// Interp value using coefficient for env info
	float_t			interp_xy(size_t index, 
								float_t x, float_t y,
								const matrix_f& coef);
public:
	/// Set the 2D axis and the curve
	void            set_data(axis_type* axis, curve_type* curve,
								const matrix_n& unknown_mat)
	{
		_axis=axis; _curve=curve; _umat=unknown_mat;
	}
	matrix_n&		get_umat() { return _umat; }
	const size_t	ghosts_count() { return _ghosts.size(); }
	const size_t	bnodes_count() { return _bnodes.size(); }
	npos_vec_type&	ghosts_npos() { return _ghosts; }
	npos_vec_type&	bnodes_npos() { return _bnodes; }
	/// Gen ghost cell info
	void			gc_gen();
protected:
	pos_vec_type	_gpos;
	pos_vec_type	_spos;
	pos_vec_type	_mpos;			// mirror surface pos
	vector_f		_marc;			// mirror arc on surface
//	pos_vec_type	_mgpos;
//	vector_n		_ms_npos;		// mirror surface npos
	/// 
	void			_prepare_pos();
public:
	pos_vec_type& 	ghosts_pos() { return _gpos; }
	pos_vec_type& 	mirrors_pos() { return _mpos; }
//	pos_vec_type&	mghosts_pos() { return _mgpos; }
	pos_vec_type& 	curve_pos() { return _spos; }
	vector_f&		mirrors_arc() { return _marc; }
//	vector_n&		msurf_npos() { return _ms_npos; }	
}; 	// end of class TFJGCGen

}	// end of namespace

#include "fjapp_gc.cpp"

#endif 
