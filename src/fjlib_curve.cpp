#include "fjlib_curve.h"
#include <cmath>
#include "fjlib_vecmat_print.h"

namespace fjlib {

std::ostream& operator<<(std::ostream& out, const TFJCurve& c)
{
	using std::endl;
	out << c.dim() << endl; 
	for (size_t i=0; i<c.dim(); i++)
	{
		vprint_f prt(c.get_points().at(i),vmpfText);
		out << prt << endl;
	}
	return out;
}

std::istream& operator>>(std::istream& in, TFJCurve& c)
{
	size_t dim;
	in >> dim;
	for (size_t i=0; i<dim; i++)
	{
		vprint_f prt(c.get_points().at(i));
		in >> prt;
	}
	return in;
}

void TFJCurve::calc_arcs()
{
	size_t n=pt_count();
	arcs_data.resize(n);
	if (n<1) return;	// throw "need at least one node";
	arcs_data[0]=0;

	if (n<2) return;
	for (size_t i=1; i<n; i++)
	{
		// calc the distance between i-1 to i nodes
		float_t sum=0.0;
		for (size_t j=0; j<dim(); j++)
			sum+=sqr(at(j,i)-at(j,i-1));
		arcs_data[i]=arcs_data[i-1]+std::sqrt(sum);
	}
}


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

}	// end of namespace
