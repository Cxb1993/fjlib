#ifndef FJAPP_DIFFUSION_GC
#define FJAPP_DIFFUSION_GC

#include "fjlib_curve.h"
#include "fjlib_scale.h"
#include "fjapp_vof.h"
#include "fjapp_gc.h"

// disable parallel by comment this line
#include "fjapp_kspm.h"

namespace fjlib {

/*!
// DiffGC simulator
// simulate a simple diffusion case where surfactant diffuses
// from a spherical surface into surroundings in axisysmetric
// coordinates. It has dirichit b.c. set at ghost point which
// makes the solution parture from the analytic solution
//
//	1.10.06 created, copied lots from fjapp_bubbleNS
//			call finalize() before exit the program
!*/

class TFJDiffGC {
public:
	typedef TFJScale_UniMesh		mesh_type;
	typedef TFJAxis_UniMesh			axis_type;
	typedef TFJCurve_CubicSpline	curve_type;
	typedef TFJVOFGen_Full<mesh_type,curve_type>
									vof_type;
	typedef TFJGCGen<mesh_type,curve_type>
									gc_type;
//private:
public:	// for testing
	float_t			_npl;		///<- nodes per unit length
	float_t			_nlen;		///<- geometry size
	float_t			_bradius;	///<- bubble radius;
	float_t			_Pe;		///<- bulk Pe number
	float_t			_acc;		///<- accuracy control
	float_t			_max_iter;	///<- max step loop iteration
protected:
	gc_type			_gc;		///<- gc object
	vof_type		_vof;		///<- vof object
	matrix_n		_pre_umat;	///<- unknown map
	axis_type		_axis;		///<- coordinates
	vector_f		_x,_y;		///<- curve coordinates
	curve_type		_surf;		///<- curved geometry
	matrix_f		_conc;		///<- bulk concentration
	matrix_f		_csu,_cap,_cae,_caw,_can,_cas;
private:
	void			_set_ghosts();
public:
	TFJDiffGC():_Pe(2.0),_npl(2),_nlen(5),_bradius(1),
				_acc(1e-12),_max_iter(500) {}
//	~TFJDiffGC() { _finalize_ksp(); }
	virtual
	void			initialize();
	virtual
	void			solve();
	virtual
	void			prepare(bool ghost_only=false);
	matrix_f&		get_cmat() { return _conc; }
#ifdef _Petsc_KSP
protected:
	TFJKSP_Mask		_cksp;
	int				_rank;
//	void			_finalize_ksp() { _cksp.finalize(); }
public:
//	void			finalize() { _finalize_ksp(); }
	void			initialize_ksp();
	void			set_rank(int __rank) { _rank=__rank; }
	const int		rank() const { return _rank; }
#endif

};	// end of class

}	// end of namespace



#endif
