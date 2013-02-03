#include "fjlib_solver.h"

namespace fjlib {

int TFJTriDag1D::solve()
{	// modified from Numerical Recipes
	size_t cn=(*fi).size();
	float_t bet;
    vector_f gam(cn);
	bet=(*ap)[0];
	if (bet == 0.0) throw "error 1 in TriDag";
	(*fi)[0]=(*sr)[0]/bet;
	for (size_t j=1;j<cn;j++) {
		gam[j]=(*ar)[j-1]/bet;
		bet=(*ap)[j]-(*al)[j]*gam[j];
		if (bet == 0.0) throw "error 2 in TriDag";
		(*fi)[j]=((*sr)[j]-(*al)[j]*(*fi)[j-1])/bet;
	}
	for (int j=cn-2;j>=0;j--)
		(*fi)[j] -= gam[j+1]*(*fi)[j+1];
	return 0;
}

int TFJCTriDag1D::solve()
{	// modified from Numerical Recipes
	float_t alpha,beta,gamma;
	size_t cn=(*fi).size();
	if (cn <= 2) throw "\nsystem too small from solve() in fjlib_solver.cpp\n";
	vector_f bb(cn),u(cn),z(cn);
	gamma = -(*ap)[0];
	bb[0]=(*ap)[0]-gamma;
	alpha=(*ar)[cn-1];
	beta=(*al)[0];
	bb[cn-1]=(*ap)[cn-1]-alpha*beta/gamma;
	for (size_t i=1;i<cn-1;i++) bb[i]=(*ap)[i];
	vector_f *bap=ap, *bsr=sr, *bfi=fi;
	ap=&bb;
	TFJTriDag1D::solve();
//	tridag(a,bb,c,r,x);
	u[0]=gamma;
	u[cn-1]=alpha;
	for (size_t i=1;i<cn-1;i++) u[i]=0.0;
	sr=&u; fi=&z;
	TFJTriDag1D::solve();
//	tridag(a,bb,c,u,z);
	ap=bap; sr=bsr; fi=bfi;
	float_t fact=((*fi)[0]+beta*(*fi)[cn-1]/gamma)/
		(1.0+z[0]+beta*z[cn-1]/gamma);
	for (size_t i=0;i<cn;i++) (*fi)[i] -= fact*z[i];
	return 0;
}

}	// end of namespace
