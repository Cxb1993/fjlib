#ifndef FJAPP_SPHERE_DIFFUSION_H
#define FJAPP_SPHERE_DIFFUSION_H

#include "fjlib_transient.h"
//#include "fjapp_SurfTsnDyn.h"
#include "fjapp_gc.h"
#include "fjlib_solver.h"
#include "fjapp_vof_staggered.h"

// disable parallel by comment this line
#include "fjapp_kspm.h"

namespace fjlib {

struct TFJSphereAdsorpParam {
/*	// obsolete
	float_t	Pe_b;			///<- Bulk Peclet number, based on desorption velocity scale
	float_t Pe_s;			///<- Surface Peclet number
	float_t k;				///<- surfactant isoterm
	float_t lamda;			///<- diffusion over desorption time scale
	float_t x_eq;			///<- Initial bulk corresponding coverage
*/
	float_t hob;			///<- h/b
	float_t Dsb;			///<- Ds/Db
	float_t lamda;			///<- tD/(1/alpha)
	float_t	x_eq;			///<- x_eq=k/(1+k)
	
	float_t	x0;				///<- Initial surface coverage
	float_t radius;			///<- Radius of the sphere
	float_t	npl;			///<- Node per unit length
	float_t	size;			///<- Domain size

	float_t dt;				///<- Time step
	float_t tmax;			///<- Max time 
	float_t alpha;			///<- damping factor
	float_t conc_acc;		
	float_t iter_acc;
	float_t	max_iters;		///<- Maximium iterations allowed

	int		debug;			///<- Debug mode
};

/*!
// SphereAdsorp simulator
// simulates a diffusion ad(de)sorption from far field to 
// a spherical interface.
// This is transient problem(static shape) and diffusion and
// ad(de)sorption equation is highly coupled at interface and
// the coupling is non-linear.
//
// ghost cell method is adapted here to demonstrated how to 
// use regular mesh to deal with irregular domain
//
// 1.20.06 created, copied most from fjapp_diffGC and GCC
// 			fjapp_bubbleNSG, fjapp_bubbleNS
// 1.30		interp_cghost() is replace by a mirror interpolation
// 2.21		modify to include staggered mesh support
!*/

class TFJSphereAdsorp: public TFJTransient {
public:
	typedef TFJSphereAdsorpParam			params_type;
	typedef TFJCurve_CubicSpline			curve_type;
	typedef TFJCurve_Line					lcurve_type;
/*
	typedef TFJScale_UniMesh				mesh_type;
	typedef TFJAxis_UniMesh					axis_type;
	typedef TFJVOFGen_Full<mesh_type,curve_type>
											vof_type;
*/
// staggered version 
	typedef TFJScale_Staggered_UniMesh		mesh_type;
	typedef TFJAxis_Staggered_UniMesh		axis_type;
	typedef TFJVOFGen_Staggered<mesh_type,curve_type>
											vof_type;				
											
	typedef TFJGCGen<mesh_type,curve_type>	gc_type;
	typedef gc_type::node_mat_type			env_type;
	typedef TFJCubicSpline					interp_type;
protected:
	params_type		_params;		///<- parameters
//	TFJSurfTsnParams
//					_surfc_params;	///<- surfactant params recasted
//	TFJSurfTsnEquil	_surfEq;		///<- static surfactant solver
//	TFJSurfTsnDyn	_surfDyn;		///<- dynamic surfactant solver
	
	vof_type		_vof;			///<- vof calculator
	axis_type		_axis;			///<- axis object
	gc_type			_gc;			///<- ghost cell handler
	matrix_n		_pre_umat;		///<- umat from vof
	vector_f		_x,_y;			///<- coordinates of curve
	curve_type		_surf;			///<- surface object
	lcurve_type		_lsurf;			///<- line version of curve
	matrix_f		_conc,_pconc;	///<- bulk concentration storage
	matrix_f		_csu,_cap,_cae,	///<- bulk solver coefficients
					_caw,_can,_cas;
	vector_f		_cghost;		///<- ghost bulk concentration
	
	vector_f		_csub;			///<- sublayer concentration
	vector_f		_cfluxsub;		///<- sublayer normal conc flux
//	vector_f		_cfluxsub_x,_cfluxsub_y;
									///<- sublayer flux at r,z direction
//	matrix_f		_cgrad_x,_cgrad_y;
									///<- diverent of conc
	env_type		_genv,			///<- mirror ghost env info
					_senv;			///<- surface node env info
	matrix_f		_gcoef,			///<- ghost cell interp coef
					_scoef;			///<- surface node interp coef
	
	vector_f		_gama,_pgama;	///<- surface concentration
	vector_f		_sigma;			///<- surface tension
	float_t			_max_gama;		///<- maximium surface tension
	vector_f		_gcoefA,_gcoefB,///<- surface solver coefficients
					_gcoefC,_gcoefR;
	float_t			_tracked_mass;	///<- total surface mass tracked

//	vector_f		_mgama;			///<- mirror gama conc
	vector_f		_mcsub;			///<- mirror sublayer conc
	vector_f		_mcfluxsub;
	gc_type::pos_vec_type
					_mgpos;			///<- mg node pos
	vector_f		_mgconc;		///<- mg node conc
	interp_type		_sc_interp;		///<- surf conc interp

	size_t			_iters;			///<- current iteration
	float_t			_tc;			///<- current time

	TFJTriDag1D		_gama_solver;	///<- tridiagonal solver
	vector_f		_xsp,_ysp;		///<- surface derivative
	gc_type::pos_vec_type
					_cn1;			///<- normal node
	vector_f		_cn1v;			///<- normal cn1 value
	
	vector_f		_mxsp,_mysp;	///<- mirror direvitive
	gc_type::pos_vec_type
					_mn1;			///<- mirror normal node
	vector_f		_mn1v;			///<- mirror nornal value

	float_t			_get_cghost_err(const vector_f& vb);
	float_t			_get_csub_err(const vector_f& vb);
private:
	vector_n		_g_src,_s_src;
//	vector_f		_tcghost;
//	matrix_f		_tconc;
	vector_f		_tgama;
	void			_set_ghosts();
	float_t			_dist(float_t x1, float_t y1,
						float_t x2, float_t y2);
//	void			_get_ghosts(const matrix_f& m, vector_f& v);
//	size_t 			_index(size_t i, int di, bool y);
//	void			_set_grad(size_t i, size_t j, bool check=false);
//	void			_restrict_vec(vector_f& v);
//	void 			_update_surf_bc(vector_f& v);
	void			_dampen_vec(const vector_f& ov,
						vector_f& nv,float_t alpha);
	vector_n		_inner_iters;
protected:
	void			prepare_eqn_gama(bool rhs_only=false);
	void			prepare_eqn_conc(bool ghost_only=false);
//	void			update_grad_conc(bool tag2_only=false);
	void			solve_gama();
	void			update_sublayer_conc();
	void			initialize_interp();
	void			interp_mirror_csub();
	void			interp_cghost();
	void			update_ghost_grad();
//	void			interp_cgrad(const matrix_f& m,vector_f& v);
	void			interp_sublayer_cflux();
//	void			update_fake_bc(const matrix_f& m);
	void			solve_conc();
public:
	void			initialize();
	void			reset();
	void			step();
	bool			is_end() { return _tc>=_params.tmax; }
	params_type&	get_params() { return _params; }
	vector_f&		gama() { return _gama; }
	matrix_f&		conc_bulk() { return _conc; }
	matrix_n&		unknown_map() { return _gc.get_umat(); }
//	vector_f&		conc_ghosts() { return _cghosts; }
	vector_f&		conc_sublayer() { return _csub; }
	vector_f&		cflux_sublayer() { return _cfluxsub; }
	size_t			oiter_count() { return _inner_iters.size(); }
	vector_n		iiter_counts() { return _inner_iters; }
	float_t			tc() { return _tc; }
protected:	
	inline size_t	xn() { return _axis.x().seg_count(); }
	inline size_t	yn() { return _axis.y().seg_count(); }
	inline size_t	sn() { return _surf.pt_count(); }
	inline size_t 	gn() { return _gc.ghosts_count(); }
	inline float_t	dx() { return _axis.x().seg_len(0); }
	inline float_t 	dy() { return _axis.y().seg_len(0); }
#ifdef _Petsc_KSP
protected:
	TFJKSP_Mask		_cksp;
public:
	void			initialize_ksp();
#endif
};	// end of class



}	// end of namespace

#endif
