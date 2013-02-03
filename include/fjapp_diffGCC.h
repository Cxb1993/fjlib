#ifndef FJAPP_DIFFUSION_GCC_H
#define FJAPP_DIFFUSION_GCC_H

#include "fjapp_diffGC.h"

namespace fjlib {

/*!
// DiffGC_Coupled
// use ghost cell method to completly solve a simple case
// where surfactant diffuses from a spherical surface into
// surroundings with dirichit b.c set at surface
//
// This name is probably incorrect, 'cause it doesn't couple
// with anything else
//
// 1.18.06 created, start from fjapp_diffGC
//
!*/

class TFJDiffGC_Coupled: public TFJDiffGC
{
/*
public:
	typedef TFJDiffGC::gc_type::node_mat_type	node_mat_type;
	typedef TFJDiffGC::gc_type::npos_vec_type	npos_vec_type;
	typedef TFJDiffGC::gc_type::pos_vec_type	pos_vec_type;
*/
public:
	typedef gc_type::node_mat_type		 env_type;
protected:
	vector_f		_csub;		// sublayer concentration
	env_type		_genv;		// surf env info
	matrix_f		_gcoef;		// fitting surface env coef

	size_t			_iters;		// iteration count;
	size_t			_ii,_jj;	// maximum error
	float_t			_tol;		// tolerance
protected:
	/// Initialize interpolation
	void			initialize_interp();
	/// Prepare coefficient for interpolation
	void 			prepare_interp();
	/// In(ex)tropolate bulk values(with surf value) to ghost nodes
	void			interp_gnodes();
	/// check error
	float_t			_check_err(const matrix_f& p, const matrix_f& q);
public:
	TFJDiffGC_Coupled(): TFJDiffGC(),_tol(1e-2) {}
	void			set_tolerance(float_t tol)
	{ _tol=tol; }
	void			initialize();
	void			solve();
	void			solve_Dirichlet();
	size_t			get_iters() { return _iters; }
};	// end of class 

}	// end of namespace

#endif
