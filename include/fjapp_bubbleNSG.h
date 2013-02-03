#ifndef TFJLIBBubbleNSGH
#define TFJLIBBubbleNSGH

#include "fjapp_bubbleNS.h"
#include "fjapp_SurfTsnDyn.h"
#include "fjlib_solver.h"

namespace fjlib {

/*!
// BubbleNSG parameters \n
// 6.24 created
*/
struct TFJBubbleNSG_Params {
	/* physical parameter */
	float_t	Pe;				// surface Peclet number
	float_t Bi;				// Biot number
	TFJSurfTsnParams
			Surfc;			// surfactant parameters
	float_t	gacc;			///<- accuracy control of gama
	float_t	x0;				// initial(eq) surfactant coverage
};

#ifndef _Surfactant
#define _Surfactant
#endif

/*!
// TFJBubbleNSG_Solver 
//	includes surfactant mass convection-diffusion eqn
// 
// 6.24 created
// 7.27 second order in time, added
// 10.23 add get_absorbed_mass()
//		move _tracked_mass from private to protected
// 10.28 add #define _Surfactant;
// 2.22.06 add an option to prepare_eqn_gama(rhs_only) for
//			compatibility with TFJBubbleNSGD_Solver
//		   change prepare_eqn_gama and solve_gama to virtual 
//			move solve_surface_mass() to step()
//			delete solve_interface() derivative
// 2.28.06 change get_total_mass to get_curve_mass()
!*/
class TFJBubbleNSG_Solver: public TFJBubbleNS_Solver {
public:
	typedef TFJBubbleNSG_Params			params_surf_type;
protected:
	params_surf_type
					_gpms;
	TFJSurfTsnEquil _surfEq;
	TFJSurfTsnDyn	_surfDyn;
	vector_f		_gama,_pgama;				// surface concentration
	vector_f		_sigma;				// surface tension
	float_t			_max_gama;
	/// Use linear interpolation for surface concentration
	TFJInterp_Line	_gInterp,_sigInterp;	
	void			remesh_gama();
	void			initialize_interface();
	void			reset_interface();
	void			remesh_interface();
//	void			solve_interface();
	void			step();
	virtual void	solve_surface_mass();
	virtual void	update_surface_tension();

	/// Solver storage for eqn gama
	vector_f		_coefA,_coefB,_coefC,_coefR;
	/// Prepare eqn coefficient for surface mass 
	virtual void	prepare_eqn_gama(bool rhs_only=false);
	/// Prepare rhs coeff depends on adsorption limit or mix
	virtual void	prepare_eqn_gama_rhs();
	void			prepare_eqn_gama2(); // Crank-Nicolson
	virtual void	solve_gama();
	/// Evaluate rhs of the surface mass eqn
	virtual float_t	get_absorbed_mass();
	virtual void	update_tracked_mass() 
	{	_tracked_mass+=get_absorbed_mass(); }
	inline
	float_t			coverage() { return _surfEq.get_Gama(); }
	inline
	float_t			c0oa() { return _surfEq.get_C(); }

	float_t			_tracked_mass;
private:
	/// scaled gama
	float_t			gama_HatD2N(float_t val,
								float_t upstream,
								float_t downstream);
	float_t			gama_HatN2D(float_t val,
								float_t upstream,
								float_t downstream);
	/// interpolate gama using higher order scheme
	float_t			gama_H2F(float_t val);
	TFJTriDag1D		_gama_solver;
public:
	params_surf_type&
					get_surf_params() { return _gpms; }
	vector_f&		get_gama_vec() { return _gama; }
	vector_f&		get_sigma_vec() { return _sigma; }
	float_t			get_tracked_mass() { return _tracked_mass; }
	virtual float_t	get_curve_mass();
	virtual float_t	get_curve_mass_err() 
	{ return (get_curve_mass()-_tracked_mass)/_tracked_mass; }
};

}	// end of namespace

#endif
