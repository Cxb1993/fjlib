#include "fjlib_integration.h"

namespace fjlib {

float_t TFJIntegrate_Trapzd::integrate(float_t x0, 
									   float_t x1, int n)
{
	float_t h=(x1-x0)/(double)n;
	float_t result=value_at(&x0)/2.0;
	float_t x=x0;
	for (size_t i=1; i<n; i++)
	{
		x+=h;
		result+=value_at(&x);
	}
	result+=value_at(&x1)/2.0;
	result*=h;
	return result;
}

float_t TFJIntegrate_GaussQuad::integrate(float_t x0, 
									   float_t x1, int o)
{
	gauss.set(x0,x1,o);
	vector_f absc=gauss.get_abscissas();	
	for (size_t i=0; i<absc.size(); i++)
		absc(i)=value_at(&absc(i));
	return gauss.apply(absc);
}


} // end of namespace