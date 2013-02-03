#include "fjapp_SurfTsnEQ.h"

namespace fjlib {

float_t TFJSurfTsnEquil::calc_Tension(const float_t x)
{ 
	return 1.0+m_Params.El*
			(std::log(1.0-x)-
			0.5*m_Params.Vir*x*x);
}

float_t TFJSurfTsnEquil::calc_BulkConc(const float_t x)
{
	return x*std::exp(m_Params.Vir*x)/(1-x);
}

void TFJSurfTsnEquil::m_NewtonStateEqn(
		void *obj, float_t *x, vector_f* v)
{
	TFJSurfTsnEquil *a=(TFJSurfTsnEquil*)obj;
	TFJSurfTsnParams& p=a->get_params();

	float_t k=p.Vir,c=a->get_C();
	float_t xp=*x;
	float_t tmp=std::exp(k*xp);
			
	(*v)[0]=xp*tmp+c*(xp-1);
	(*v)[1]=tmp*(1+k*xp)+c;
}

void TFJSurfTsnEquil::m_NewtonTensionEqn(
		void *obj, float_t *x, vector_f* v)
{
	TFJSurfTsnEquil *a=(TFJSurfTsnEquil*)obj;
	TFJSurfTsnParams& p=a->get_params();

	float_t k=p.Vir, ela=p.El, xp=*x;
	
	(*v)[0]=1.0+ela*(std::log(1.0-xp)-0.5*k*xp*xp);
	(*v)[1]=-ela*(1.0/(1.0-xp)+k*xp);
}

void TFJSurfTsnEquil::solve_Gama()
{
	root.set_func(m_StateProc);
	root.set_range(0.0,1.0);
//	root.set_eps(1e-5);
	root.solve();
	if (!root.empty())
		_x=root.root(0);
	else
		throw "no root in solve_Gama() at fjapp_SurfTsnEQ.cpp";
}

float_t TFJSurfTsnEquil::solve_MaxGama()
{
	root.set_func(m_TensionProc);
	root.set_range(0.0,0.9999999);
	root.solve();
	if (!root.empty())
		return root.root(0);
	else
		throw "no root in solve_MaxGama() at fjapp_SurfTsnEQ.cpp";
}


}	// end of namespace