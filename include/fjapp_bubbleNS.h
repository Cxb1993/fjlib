#ifndef TFJLIBBubbleNSH
#define TFJLIBBubbleNSH

#include "fjlib_transient.h"
#include "fjlib_curve.h"
#include "fjlib_scale.h"
#include "fjlib_smoothing.h"
#include "fjapp_vof_staggered.h"

// disable parallel by comment this line
#include "fjapp_ksp.h"

namespace fjlib {

/*!
// BubbleNS parameters \n
// including physical and simulation parameters
//	
*/
struct TFJBubbleNS_Params	{
	/* physical parameter */
	float_t Re;				///<- Renolds number
	float_t Ca;				///<- Capillory number
	float_t Bo;				///<- Bond number
	float_t	Viscosity_Ratio;///<- viscous ratio of internal over external
	float_t Density_Ratio;	///<- density ratio of internal over external
	float_t Pn;				///<- pressure scaling, 1.0
	float_t	IPos;			///<- needle position, 1.0
	float_t	Uc;				///<- inject speed, 1.0

	/* simulation parameter */
	float_t	xnpl,ynpl;		///<- nodes per unit len in x&y mesh
	float_t	snpl;			///<- nodes per unit len on the surface
	bool	flat;			///<- start from flat or curved surface
	float_t xmax,ymax;		///<- 2d domain size
	float_t	uvacc,pacc;		///<- accuracy control of each variable
	float_t acc,acc_abs;	///<- accuracy of the whole solver		
	float_t	dt,dt_save;		///<- simulation time step
	float_t	tmax;			///<- max simulation time step
	float_t beta;			///<- pressure splitting coefficient
	size_t	max_iter1;		///<- 1st loop iteration max steps for NS
	size_t	max_iter2;		///<- 2nd loop iteration max steps for NS
	float_t	t0;				///<- start time step
};

/*!
// BubbleNS Interface class \n
// used to manage data and functions assorciated with the interface
//
// Cation:
//  It takees two steps to remesh: 1. call remesh_??? to generate
//  a new mesh; 2. after done all the remesh work related to other
//  properties, call finalize_mesh()
//
//
//	6.24 Add finalize_mesh() to make it possible that
//		other properties associated with nodes can be splined
//		simultaneously
//	6.27 Add calc_curvature() 
//	10.20 need to change to be able to identify the last 
//		interface and apply bc before spline and u,v bc 
//		add exposure of internal _x _y 
//		add type_id and set_type() to reflect the change due to
//			different part of interface
//			type_id=0 main drop, bottom one
//			type_id=1 rising drop, middle or up one
*/
class TFJBubbleNS_Interface {
public:
	typedef TFJCurve_CubicSpline	curve_type;
	typedef TFJAxis_Staggered_Mesh	axis_type;
protected:
	curve_type		_curve;		// curve object
	float_t			npl;		// nodes per unit length
	vector_f		_x,_y;		// nodes pos storage
	vector_f		_ns;		// arc vector for new mesh
	bool			flat;
	void			initialize();
public:
	inline
	float_t			get_npl() { return npl; }
	inline
	bool			get_flat() { return flat; }

	curve_type&		curve() { return _curve; }
	const
	curve_type&		curve() const { return _curve; }
	/// initialize curve
	virtual
	void			reset();
	/// spline the whole curve and its related properties
	virtual
	void			spline();	
	virtual
	void			apply_bc();
	virtual
	void			calc_derivs();
	virtual
	void			arc_deriv(const vector_f& v, size_t ord, 
							vector_f& ov);
	/// remesh is used to make a new mesh
	virtual
	bool			remesh();
	virtual
	void			remesh_uniform();
	virtual
	void			remesh(const axis_type& axis);
	vector_f&		new_mesh_arc() { return _ns; }
	/// finalize_mesh copies new mesh over
	virtual
	void			finalize_mesh() { _curve.remesh(_ns); }
public:
	virtual
	void			set(float_t nodes_per_len, bool is_flat=true)
	{ npl=nodes_per_len; flat=is_flat; initialize(); }
	TFJBubbleNS_Interface(): npl(1.0), flat(true)
	{ initialize(); } 
protected:		// make them public for convenience
	vector_f		_xp1,_xp2;	// derivative of x(s)
	vector_f		_yp1,_yp2;	// derivative of y(s)
	vector_f		_arc_crt;	// arc length correction
public:
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
	vector_f&		arc_crt() { return _arc_crt; }
	virtual 
	void			calc_curvature(vector_f &cvt);
public:
	inline
	vector_f&		x() { return _x; }
	inline
	vector_f&		y() { return _y; }
	inline
	const vector_f&	x() const { return _x; }
	inline
	const vector_f&	y() const { return _y; }
protected:
	size_t			type_id;	
public:
	inline
	void			set_type(size_t type) 
	{ 
		type_id=type;
		apply_spline_bc(); 
	}
	inline
	size_t			get_type() { return type_id; }
	void			apply_spline_bc();
};

/*!
//	5.29 add helper functions to manuplate axis
//	5.30 find a bug in diffusion term in backup_uv, add /2.0 
//		but no much impact so far
//	6.25 5.30 bug is not be a bug, over-corrected, recovered
//	7.1 find a bug change to _pu instead _u
		switch to second order scheme in time
//  10.19 add upate_surf_vel_bc
//	1.9.06	move parallel solver support here in case of not
//			need to mantain multiple _ksp.cpp files
//	2.24.06 find a bug for non-flat start up, weird!!
//			find out it's due to interpolation.cpp 1e-16
// Caution:
//	_uap*Up+_uae*Ue+_uaw*Uw+_uan*Un+_uas*Us+Usu=0
//	the above applies to v and p eqnation
//
//	the interface is an object wraps a curve object
//	the axis object is a staggered helper object, use axis(sptO) etc. 
//	
!*/
class TFJBubbleNS_Solver: public TFJTransient {
public: 
	typedef TFJScale_Staggered_UniMesh	mesh_type;
	typedef TFJAxis_Staggered_UniMesh	axis_type;
//	typedef TFJScale_Staggered_Mesh		mesh_type;
//	typedef TFJAxis_Staggered_Mesh		axis_type;
	typedef TFJBubbleNS_Interface		surface_type;
	typedef surface_type::curve_type	curve_type;
	typedef TFJVOFGen_Staggered_ST<mesh_type,curve_type>
										vof_type;
	typedef TFJStaggeredAxis_Helper<mesh_type>
										staggered_axis_type;
	typedef TFJStaggeredPosType			staggered_type;
	typedef TFJBubbleNS_Params			params_type;
public:
	void			initialize();
	void			reset();
	void			step();
	void			after_step();
	bool			is_end() { return (_tc>_pms.tmax); }
protected:
	float_t			_tc;
	/// all the parameters
	params_type		_pms;
	// interface object
	surface_type	_surf;
	///	interface velocity
	vector_f		_us,_vs;	// velocity of nodes
	vector_f		_un,_ut;	// normal and tangiential velocity of nodes
	/// curvature, reserved for surface mass balance
	vector_f		_cvt;
	/// vof object for staggered mesh
	vof_type		_vof;
	// axis for the staggered mesh
	staggered_axis_type
					_axis;
//	axis_type		_axis,_uaxis,_vaxis,_paxis;
	// unkown variables
	matrix_f		_u,_v,_p;
	// unkown from last time step;
	matrix_f		_pu,_pv,_pp;
	// solver storage for eqn u,v,p
	matrix_f		_usu,_uap,_uae,_uaw,_uan,_uas;
	matrix_f		_vsu,_vap,_vae,_vaw,_van,_vas;
	matrix_f		_psu,_pap,_pae,_paw,_pan,_pas;
	// convective and surface tension terms at n and n-1 time step
	matrix_f		_ucv,_upcv,_uppcv,_vcv,_vpcv,_vppcv;
	matrix_f		_ust,_upst,_uppst,_vst,_vpst,_vppst;
	// diffusion terms at n time step
	matrix_f		_updf,_vpdf;
	///	vof values for staggered mesh
	matrix_f		_uvof,_vvof,_pvof,_nvof;
	/// residue matrix
	matrix_f		_res;
	/// iteration count
	size_t			niter,piter;
protected:
	///	solve navier stokes eqn u,v,p
	virtual void	solve_ns();
	/// solve the new interface and its related
	virtual void	solve_interface();
	/// backup for 1st or 2nd order step scheme
	void			prepare_step_scheme(bool second);
	/// prepare data for 2nd order step scheme u
	void			backup_data_u();
	/// prepare data for 2nd order step scheme v
	void			backup_data_v();

	/// prepare the eqnations for u
    void			prepare_eqn_u(bool p_only=false);
	/// prepare the eqnations for v
	void			prepare_eqn_v(bool p_only=false);
	/// prepare the eqnations for p
	void			prepare_eqn_p(bool uv_only=false);

	/// correct the bc for equations u,v,p
	virtual void	prepare_bc(staggered_type st) {}

	/// correct uv to uv star
	void			correct_uv(float_t sign);
	/// update u to boundary nodes
	void			update_bc_u();
	/// update v to boundary nodes
	void			update_bc_v();
	/// interpolate velocity for the interface
	void			get_interface_uv();
	/// initialize interface
	virtual void	initialize_interface()
	{ 	_surf.set(_pms.snpl,_pms.flat); }
	/// reset interface
	virtual void	reset_interface();
	/// spline interface
	virtual void	spline_interface() { _surf.spline(); }
	/// remesh interface
	virtual void	remesh_interface();
	///	move interface to new location
	void			move_interface();
	///	prepare interface for the new location
	void			prepare_interface();
	/// calculate the vof on the mesh
	void			update_vof_st();
	/// calculate the surface force for u,v eqn
//	void			update_surface_force();
	/// CGSTAB solver
	bool			CGSTAB(size_t xc, size_t yc,
							float_t acc, size_t niter,
							matrix_f& x, matrix_f& ap,
							matrix_f& ae, matrix_f& aw,
							matrix_f& an, matrix_f& as,
							matrix_f& b, matrix_f& res);
							
protected:		/// helper function
	///	helper, caluclate the max absolute value of the matrix elements
	float_t			mat_max_abs(matrix_f& ma);
	void			fill_vof_mat(staggered_type stag_loc,
									matrix_f &mat);
	void			fill_st_mat(staggered_type stag_loc,
									bool fill_ust, matrix_f &mat);
	inline
	float_t			den(float_t vof_v)
	{ return 1.0-(1.0-_pms.Density_Ratio)*vof_v; }
	float_t			vis(float_t vof_v)
	{ return 1.0-(1.0-_pms.Viscosity_Ratio)*vof_v; }
public:
	const
	staggered_axis_type&
					get_axis() const { return _axis; }
	params_type&	get_params() { return _pms; }
	surface_type&	get_surface() { return _surf; }
	float_t			tc() const { return _tc; }

	matrix_f&		get_vof_mat() { return _pvof; }
	matrix_f&		get_u_mat() { return _u; }
	matrix_f&		get_v_mat() { return _v; }
	matrix_f&		get_p_mat() { return _p; }
	surface_type::curve_type&
					get_curve() { return _surf.curve(); }
	vector_f&		get_un_vec() { return _un; }
	vector_f&		get_ut_vec() { return _ut; }
	vector_f&		get_us_vec() { return _us; }
	vector_f&		get_vs_vec() { return _vs; }
	vector_f&		get_cvt_vec() { return _cvt; }

	const matrix_f&	get_vof_mat() const { return _pvof; }
	const matrix_f&	get_u_mat() const { return _u; }
	const matrix_f&	get_v_mat() const { return _v; }
	const matrix_f&	get_p_mat() const { return _p; }
	const surface_type::curve_type&
					get_curve() const { return _surf.curve(); }

	const size_t	get_iter1() const { return niter; }
	const size_t	get_iter2() const { return piter; }

	float_t			get_tracked_vol();
	float_t			get_total_vol() { return M_PI*_tc; }
	float_t			get_vol_err() 
	{ 
		float_t vol=get_total_vol();
		return (get_tracked_vol()-vol)/vol;
	}
public:
	float_t			save_time;		// temperory
public:
	// helper function for fast access with less coding
	inline const
	float_t			dx(int m) const { return _axis.dx(m); }
	inline const
	float_t			dy(int n) const { return _axis.dy(n); }
	inline const
	float_t			dx2(int m) const { return _axis.dx(m)/2.0; }
	inline const
	float_t			dy2(int n) const { return _axis.dy(n)/2.0; }
	inline const
	float_t			xp(int m) const { return _axis.xp(m); }
	inline const
	float_t			yp(int n) const { return _axis.yp(n); }
	inline const 
	size_t			xc() const { return _axis.xc(); }
	inline const 
	size_t			yc() const { return _axis.yc(); }
	inline const 
	size_t			sc() { return _surf.curve().pt_count(); }

	// more helper function with restriction
	inline const	
	float_t			xe(int m) const { return xp(m)+dx2(m); }
	inline const	
	float_t			xw(int m) const { return xp(m)-dx2(m); }
	inline const	
	float_t			yn(int n) const { return yp(n)+dy2(n); }
	inline const	
	float_t			ys(int n) const { return yp(n)-dy2(n); }
	inline const
	float_t			dxe(int m) const { return (dx(m)+dx(m+1))/2.0; }
	inline const
	float_t			dxw(int m) const { return (dx(m)+dx(m-1))/2.0; }
	inline const
	float_t			dyn(int n) const { return (dy(n)+dy(n+1))/2.0; }
	inline const
	float_t			dys(int n) const { return (dy(n)+dy(n-1))/2.0; }

public:
	virtual void	apply_interface_uv_bc();

#ifdef _Petsc_KSP
protected:
	TFJKSP			_uksp,_vksp,_pksp;
//	void			solve_ns();
	int				_rank;
public:
	void			initialize_ksp();
	void			set_rank(int __rank) { _rank=__rank; }
	const int		rank() const { return _rank; }
#endif

};


}	// end of namespace

#endif
