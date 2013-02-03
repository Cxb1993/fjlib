//---------------------------------------------------------------------------
// Base class for CFD solver

#ifndef FJLIBSolverH
#define FJLIBSolverH

#include "fjlib.h"

namespace fjlib {

/*!
// Base class for solver
//
// 1.24 revised, 2D removed
*/

class TFJSolver {
public:
	/// solve and returns the error code		
	virtual 
	int			solve()=0;
};

/*!
// Tridiagonal solver
//
// AP*FI(n)+AL*FI(n-1)+AR*FI(n+1)=SR
*/
class TFJTriDag1D: public TFJSolver {
protected:
    vector_f	*ap,*al,*ar,*sr,*fi;
public:
	void		set_data(vector_f *_ap, vector_f *_al,
						vector_f *_ar, vector_f *_sr,
						vector_f *_fi)
	{ ap=_ap; al=_al; ar=_ar; sr=_sr; fi=_fi; }
	int			solve();
};

/*!
// Cyclic Tridiagonal solver
//
// Same as tridiagonal solver except the last node 
// is connected with the first node
*/
class TFJCTriDag1D: public TFJTriDag1D {
// Cyclic Tridiagonal Systems
public:
	int			solve();
};

}	// End Of Namespace

#endif

