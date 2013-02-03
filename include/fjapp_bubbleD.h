#ifndef FJAPP_BUBBLE_DIFFUSION_H
#define FJAPP_BUBBLE_DIFFUSION_H

#include "fjlib_transient.h"
#include "fjapp_gc.h"
#include "fjapp_interface.h"
#include "fjapp_vof.h"

#include "fjapp_kspm.h"

namespace fjlib {

struct TFJBubbleD_Param {
	float_t Pe;				///<- Peclet number
	float_t	npl;			///<- Node per unit length
	float_t snpl;			///<- Surface node per unit length
	float_t	xmax,ymax;		///<- Domain size
	float_t t0,dt,tmax;		///<- Time step
	float_t conc_acc;		///<- Conc eqn accuracy
};

/*!
// Diffusion Problem in NonRegular Geometry
// regular mesh, not staggered, 'cause only conc needs to be solved
//
// 3.29.06 created, copied most from fjapp_SphereAdsorp 
//			fjapp_bubbleNS, and fjapp_bubbleNSGD
//			to suit as a simple protocol class
!*/

class TFJBubbleD_Solver: public TFJTransient {
public:
	typedef TFJBubbleD_Param				params_type;
	typedef TFJInterface					surface_type;
	typedef TFJCurve_CubicSpline			curve_type;
	typedef TFJCurve_Line					lcurve_type;

	typedef TFJScale_UniMesh				mesh_type;
	typedef TFJAxis_UniMesh					axis_type;
	typedef TFJVOFGen_Full<mesh_type,curve_type>
											vof_type;
											
	typedef TFJGCGen<mesh_type,curve_type>	gc_type;
	typedef gc_type::node_mat_type			env_type;
protected:
	float_t			_tc;			///<- current time
protected:
	params_type		_pms;			///<- parameters
	
	vector_f		_x,_y;			///<- coordinates of curve
	vector_f		_su,_sv;		///<- surface velocity
	surface_type	_surf;			///<- surface object
	axis_type		_axis;			///<- axis object
	vof_type		_vof;			///<- vof calculator
	matrix_f		_vof_mat;		///<- vof value
	gc_type			_gc;			///<- ghost cell handler
	matrix_n		_pre_umat;		///<- umat from vof
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
	vector_f		_mcfluxsub;		///<- mirror sublayer flux
	TFJInterp_Line	_sc_interp;		///<- interpolator mirror csub

	gc_type::pos_vec_type
					_cn1;			///<- normal node
	vector_f		_cn1v;			///<- normal cn1 value
	
	vector_f		_mxsp,_mysp;	///<- mirror direvitive
	gc_type::pos_vec_type
					_mn1;			///<- mirror normal node
	vector_f		_mn1v;			///<- mirror nornal value
private:
	void			_set_ghosts();
	float_t			_dist(float_t x1, float_t y1,
						float_t x2, float_t y2);
	void			_dampen_vec(const vector_f& ov,
						vector_f& nv,float_t alpha);
protected:
	void 			spline_interface() { _surf.spline(); }
	void			remesh_interface();
	void			move_interface();
	void			update_interface_velocity();
	void			prepare_eqn_conc(bool ghost_only=false);
	void			interp_mirror_csub();
	void			interp_cghost();
	void			interp_sublayer_cflux();
	void			solve_eqn_conc();
	
	void 			initialize_interface();
	void			reset_interface();
	void			solve_interface();
	void			initialize_mesh();
	void			initialize_vof();
	void			initialize_gc();
	void			prepare_vof();
	void 			prepare_gc();
	void			solve_bulk_mass();
public:
	void			initialize();
	void			reset();
	void			step();
	void			after_step() { _tc+=_pms.dt; }
	bool			is_end() { return _tc>=_pms.tmax; }
	float_t			tc() { return _tc; }
	matrix_n&		get_unknown_map() { return _gc.get_umat(); }
	curve_type&		get_curve() { return *_surf.curve(); }
	params_type&	get_params() { return _pms; }
	matrix_f&		get_bulk_conc() { return _conc; }
	vector_f&		get_sublayer_cflux() { return _cfluxsub; }
	vector_f&		get_ghost_conc() { return _cghost; }
	vector_f&		get_sublayer_conc() { return _csub; }
	matrix_f&		get_vof() { return _vof_mat; }
	vector_f&		get_su() { return _su; }
	vector_f&		get_sv() { return _sv; }
protected:	
	inline size_t	xc() { return _axis.x().seg_count(); }
	inline size_t	yc() { return _axis.y().seg_count(); }
	inline size_t	sc() { return _surf.curve()->pt_count(); }
	inline size_t 	gc() { return _gc.ghosts_count(); }
	inline float_t	dx() { return _axis.x().seg_len(0); }
	inline float_t 	dy() { return _axis.y().seg_len(0); }
#ifdef _Petsc_KSP
protected:
	TFJKSP_Mask		_cksp;
public:
	void			reset_cksp();
	void			initialize_cksp();
#endif
};	// end of class



}	// end of namespace

#endif
