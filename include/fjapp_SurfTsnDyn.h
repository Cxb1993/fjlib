#ifndef FJLib_Surface_Tension_Dynamics_H
#define FJLib_Surface_Tension_Dynamics_H
/*
	Surfactant dynamic class,
	surf conc scaled with Gama0,eq
	bulk conc scaled with C0,eq
	surf tension scaled with Sigma0,eq

	flux=(Coa*Cs*(1/x0-Gama)-exp(K*Gama*x0)*Gama)

	the flux is scaled depending on the cases,
	it can be Bi or Lamda. You have to know it,
	before plug into the calculation.

	Cs=const simulates ad(de)sorption calculation

	It has to be coupled with diffusion process to
	be able to solve dynamic surface tension revolution
	with iteration theme.
*/
#include "fjapp_SurfTsnEQ.h"

namespace fjlib {

class TFJSurfTsnDyn {
public:
	typedef	TFJSurfTsnEquil SurfTsnEq_Type;
private:
	// Surfactant Euqilibrium class	
	SurfTsnEq_Type*	m_pEq0;
	SurfTsnEq_Type	m_Eq;	// used for intermidiete calculation
	//	
//	virtual
//	void			m_NewtonAdsorpEqn(void *obj, float_t *x, vector_f* v);
public:
	void			set_Eq(SurfTsnEq_Type* eq) 
	{ 
		m_pEq0=eq; 
		m_pEq0->calc_C_Sigma();
		m_Eq.set_params(eq->get_params());
	}
	inline
	SurfTsnEq_Type*	Eq() { return m_pEq0; }
protected:
	float_t			_Gama,		// Surf conc scaled with Gama0,eq
					_C,			// Bulk conc scaled with C0,eq
					_Sigma;		// Surf tension scaled with Sigma0,eq
public:
	// helper equlibrium 
	inline
	float_t			GamaEq(bool dimensional) 
	{ return m_pEq0->get_Gama(dimensional); }
	inline
	float_t			CEq(bool dimensional) 
	{ return m_pEq0->get_C(dimensional); }
	inline
	float_t			SigmaEq(bool dimensional) 
	{ return m_pEq0->get_Sigma(dimensional); }

	// helper conversion functions
	inline
	float_t			Gama_N2D(const float_t gama)
	{ return gama*GamaEq(true); }
	inline
	float_t			Gama_D2N(const float_t gama)
	{ return gama/GamaEq(true); }
	inline
	float_t			C_N2D(const float_t c)
	{ return c*CEq(true); }
    inline
	float_t			C_D2N(const float_t c)
	{ return c/CEq(true); }
	inline
	float_t			Sigma_N2D(const float_t sigma)
	{ return sigma*SigmaEq(true); }
	inline
	float_t			Sigma_D2N(const float_t sigma)
	{ return sigma/SigmaEq(true); }
public:
	void			set_Gama(const float_t gama, bool dimensional=false)
	{
		if (!dimensional) _Gama=gama; 
		else _Gama=Gama_D2N(gama);
	}
	void			set_C(const float_t cs, bool dimensional=false)
	{
		if (!dimensional) _C=cs; 
		else _C=C_D2N(cs);
	}
	void			set_Sigma(const float_t sigma, bool dimensional=false)
	{
		if (!dimensional) _Sigma=sigma; 
		else _Sigma=Sigma_D2N(sigma);
	}
	float_t			get_Gama(bool dimensional=false)
	{
		if (!dimensional) return _Gama;
		else return Gama_N2D(_Gama); 
	}
	float_t			get_C(bool dimensional=false)
	{
		if (!dimensional) return _C;
		else return C_N2D(_C); 
	}
	float_t			get_Sigma(bool dimensional=false)
	{
		if (!dimensional) return _Sigma;
		else return Sigma_N2D(_Sigma); 
	}
public:
	void			calc_Sigma()
	{
		m_Eq.set_Gama(Gama_N2D(_Gama),true);
		m_Eq.calc_Sigma();
		set_Sigma(m_Eq.get_Sigma(true),true);
	}

/*
		inline float_t	Get_Flux()
		{ return ks()*_Cs*(1.0/x0()-_Gama)-
					std::exp(Eq0.Vir*_Gama*x0())*_Gama; }

		void			Solve_Cs(const float_t flux)
		{ _Cs=(flux+std::exp(Eq0.Vir*_Gama*x0())*_Gama)/ks()/
				(1.0/x0()-_Gama); }
        
		void			Solve_Gama(const float_t flux);

		inline float_t	Get_Max_Gama()
		{ return Eq0.Get_Max_Gama()/Eq0.Get_Gama_NonDim(); }
private:
		TFJRoot_Newton	root;
public:
		float_t			_flux;			// for root solving
*/
};

}

#endif