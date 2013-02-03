#ifndef FJLIB_BUBBLE_NSGD_SOLVER_H
#define FJLIB_BUBBLE_NSGD_SOLVER_H

#include "fjapp_bubbleNSG.h"
#include "fjapp_gc.h"
#include "fjapp_kspm.h"

#define _Bulk_Diffusion

namespace fjlib {

enum TFJRateLimitStepType {
	rlsInsoluable=0,		///	0
	rlsAdsorption,			///	1
	rlsDiffusion,			/// 2
	rlsMixed				/// 3
};

// created 2.22.06
struct TFJBubbleNSGD_Params {
	float_t hob;			///<- h/b
	float_t Peb;			///<- bulk Peclet number 
	TFJRateLimitStepType
			rate_limit_step;		
							///<- rate limiting step type
	
	float_t alpha;			///<- damping factor
	float_t conc_acc;		///<- bulk conc accuracy control 
	float_t iter_acc;		///<- inner iteration accuracy	
	float_t	max_iters;		///<- Maximium iterations allowed
};

// add diffusion support to bubble
// for prototype, check TFJSphereAdsorp 
// 2.21.06 created
// 2.27.06 add support for zero count of ghost cell
// 2.28.06 add support for no-mirror ghost node
// 3.21.06 add full controlled type support, rate_limit_step
// 3.23.06 found a bug Peb instead of Pe, fuck!
class TFJBubbleNSGD_Solver:public TFJBubbleNSG_Solver {
/*
public: 
	typedef TFJScale_Staggered_UniMesh	mesh_type;
	typedef TFJAxis_Staggered_UniMesh	axis_type;
	typedef TFJBubbleNS_Interface		surface_type;
	typedef surface_type::curve_type	curve_type;
	typedef TFJVOFGen_Staggered_ST<mesh_type,curve_type>
										vof_type;
	typedef TFJStaggeredAxis_Helper<mesh_type>
										staggered_axis_type;
	typedef TFJStaggeredPosType			staggered_type;
	typedef TFJBubbleNS_Params			params_type;
*/
public:
	typedef TFJBubbleNSGD_Params			params_diff_type;
	typedef TFJGCGen<mesh_type,curve_type>	gc_type;
	typedef gc_type::node_mat_type			env_type;
protected:
//	bool			_diffusion_switch;	///<- turn off, NSG
	params_diff_type
					_dpms;			///<- diffusion parameters

	matrix_n		_pre_umat;		///<- pre unknown matrix
	gc_type			_gc;			///<- gc object
	matrix_f		_conc,_pconc;	///<- bulk concentration storage
	matrix_f		_csu,_cap,_cae,	///<- bulk solver coefficients
					_caw,_can,_cas;
	vector_f		_cghost;		///<- ghost bulk concentration
	vector_f		_csub;			///<- sublayer concentration
	vector_f		_cfluxsub;		///<- sublayer normal conc flux
	
	env_type		_mn1env,		///<- mirror normal mn1 env info
					_cn1env;		///<- curve normal sn1 env info
	matrix_f		_mn1coef,		///<- mn1 node interp coef
					_cn1coef;		///<- cn1 node interp coef
	vector_n		_mn1src,		///<- mn1 source index
					_cn1src;		///<- cn1 source index
	
	vector_f		_mcsub;			///<- mirror sublayer conc
	vector_f		_mcfluxsub;
	TFJInterp_Line	_sc_interp;
//	TFJCubicSpline	_sc_interp;		///<- surf conc interp

	gc_type::pos_vec_type
					_cn1;			///<- normal node
	vector_f		_cn1v;			///<- normal cn1 value
	
	vector_f		_mxsp,_mysp;	///<- mirror direvitive
	gc_type::pos_vec_type
					_mn1;			///<- mirror normal node
	vector_f		_mn1v;			///<- mirror nornal value
	
	float_t			_get_cghost_err(const vector_f& vb);
	float_t			_get_csub_err(const vector_f& vb);
	vector_n		_inner_iters;
private:
	/// copy ghost value into bulk matrix
	void			_set_ghosts();
	/// calculate distance between two coord
	float_t			_dist(float_t x1, float_t y1,
						float_t x2, float_t y2);
	/// dampen new value by alpha with old value
	void			_dampen_vec(const vector_f& ov,
						vector_f& nv,float_t alpha);
	// probably obsolete
	TFJInterp_Line	_fluxInterp;	
protected:
	void			initialize_gc();
	void 			prepare_gc();
	void			prepare_eqn_gama_rhs();
//	void			prepare_eqn_gama(bool rhs_only=false);
//	void			solve_gama();
	
	void			prepare_eqn_conc(bool ghost_only=false);
	void			update_sublayer_conc();
	void			guess_sublayer_conc();
	void			initialize_env_interp();
	void			interp_mirror_csub();
	void			interp_cghost();
	void			update_ghost_grad(); // prob obsolete too
	void			interp_sublayer_cflux();
	void			solve_conc();
	/// solve surface and bulk mass balance
	virtual void 	solve_coupled_mass();
	float_t			get_absorbed_mass();
public:
//	TFJBubbleNSGD_Solver() {}
//	void			switch_diffusion_on(bool diff=true)
//	{ _diffusion_switch=diff; }
	params_diff_type&
					get_diff_params() { return _dpms; }
	inline const 
	size_t 			gc() { return _gc.ghosts_count(); }
	matrix_f&		get_conc_mat() { return _conc; }
	virtual float_t	get_bulk_mass();
	virtual float_t get_bulk_mass_err();
	vector_f&		get_sublayer_flux() { return _cfluxsub; }
	vector_f&		get_ghost_conc() { return _cghost; }
	vector_f&		get_sublayer_conc() { return _csub; }
	float_t			lamda();
public:
	void			initialize();
	void			reset();
	void			step();
#ifdef _Petsc_KSP
protected:
	TFJKSP_Mask		_cksp;
	void			reset_cksp();
	void			initialize_cksp();
#endif
};	// end of class

}	// end of namespace

#endif
