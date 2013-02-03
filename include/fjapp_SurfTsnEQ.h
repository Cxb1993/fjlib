#ifndef FJLib_Surface_Tension_Equlibrium_H
#define FJLib_Surface_Tension_Equlibrium_H

#include <cmath>
#include "fjlib_root.h"
#include <boost/bind.hpp>

namespace fjlib {


/*!
// Equalibrium Surface Tension calculation
// based on frumkin model
	surf conc scaled with GamaInf
	bulk conc scaled with a
	surf tension scaled with Sigma0
// Couple of problem can be solved, ex
//	Gama known to calculate C or C known to solve Gama
//  Gama, C known to solve for the maximum Sigma etc.
//
// Helper function conversion between non-dimensional to dimensional
// N2D, D2N
// Notes:
//	2003(4) created
//	6.20.2005 Updates this class, make it more general
!*/

// Surface model parameter based on Frumkin
struct TFJSurfTsnParams {
	float_t			GamaInf,		// Max packing surf conc
					a,				// desorption/adsorption const						
					Vir,			// Frumkin virial coef
					El,				// Surf tension elasticity
					Sigma0;			// Clean surf tension
};

class TFJSurfTsnEquil {
private:
	TFJProcIn1OutNp	m_StateProc,		// Eqation of state
					m_TensionProc;	// Eqation of tension 
	TFJRoot_Newton	root;			// Newton solver
protected:
	TFJSurfTsnParams
					m_Params;		// Surface raw properties
	// used for solving surface concentration 
	virtual
	void			m_NewtonStateEqn(void *obj, float_t *x, vector_f* v);
	// used for solving maximium surface concentration when surfac tension equals to zero
	virtual
	void			m_NewtonTensionEqn(void *obj, float_t *x, vector_f* v);
	// calcuate surface tension
	virtual
	float_t			calc_Tension(const float_t x);
	// calucate bulk concentration from eqn of state
	virtual
	float_t			calc_BulkConc(const float_t x);

	float_t			_x,				// Surf conc scaled with GamaInf
					_ks,			// Bulk conc scaled with a
					_sig;			// Surf tension scaled with Sigma0
public: 
	// conversion helper function 
	inline 
	float_t			Gama_N2D(const float_t x) 
	{ return x*m_Params.GamaInf; }
	inline
	float_t			Gama_D2N(const float_t gama) 
	{ return gama/m_Params.GamaInf; }
	inline 
	float_t			C_N2D(const float_t ks)
	{ return ks*m_Params.a; }
	inline
	float_t			C_D2N(const float_t c)
	{ return c/m_Params.a; }
	inline
	float_t			Sigma_N2D(const float_t sig)
	{ return sig*m_Params.Sigma0; }
	inline 
	float_t			Sigma_D2N(const float_t sig)
	{ return sig/m_Params.Sigma0; }
public:
	TFJSurfTsnEquil() 
	{
		// set data obj for solver
		root.set_obj(this);

		// bind to member function
		m_StateProc=boost::bind(
			&TFJSurfTsnEquil::m_NewtonStateEqn,
			this,_1,_2,_3);
		m_TensionProc=boost::bind(
			&TFJSurfTsnEquil::m_NewtonTensionEqn,
			this,_1,_2,_3);
	}
	void			set_params(TFJSurfTsnParams p)
	{ m_Params=p; }
	TFJSurfTsnParams&
					get_params() { return m_Params; }

	void			set_Gama(const float_t x, bool dimensional=false)
	{ 
		if (!dimensional) _x=x;
		else _x=Gama_N2D(x);
	}
	void			set_C(const float_t ks, bool dimensional=false)
	{
		if (!dimensional) _ks=ks;
		else _ks=C_N2D(ks);
	}
	inline
	float_t 		get_Gama(bool dimensional=false)
	{
		if (!dimensional) return _x;
		else return Gama_N2D(_x);
	}
	inline
	float_t			get_C(bool dimensional=false)
	{
		if (!dimensional) return _ks;
		else return C_N2D(_ks);
	}
	inline
	float_t			get_Sigma(bool dimensional=false)
	{
		if (!dimensional) return _sig;
		else return Sigma_N2D(_sig);
	}

public:
	// Calculate sigma from gama
	void			calc_Sigma() { _sig=calc_Tension(_x); }
	// Calculate bulk conc from gama
	void			calc_C() { _ks=calc_BulkConc(_x); }
	// Solve gama from bulk conc
	void			solve_Gama();
	// Solve max gama when sigma=0
	float_t			solve_MaxGama();
	// Calculate sigma and bulk conc from gama
	void			calc_C_Sigma() { calc_C(); calc_Sigma(); }
};

}	// end of namespace

#endif