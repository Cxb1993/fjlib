#ifndef FJLIBGAUSSQUADH
#define FJLIBGAUSSQUADH

#include "fjlib_float.h"
#include "fjlib_vector.h"

namespace fjlib {
	
class TFJGaussQuadrature {
private:
	void		gauleg(const float_t x1, const float_t x2,
						vector_f &x, vector_f &w);
protected:
	vector_f	abscissas,weights;
public:
	int			order;
	void		set(float_t low, float_t high, int o);

	vector_f&	get_abscissas() { return abscissas; }
	float_t		apply(const vector_f &v);
};

} // end of namespace

#endif
