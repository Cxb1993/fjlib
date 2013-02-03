#include "fjapp_diffGCC.h"
#include "fjlib_vecmat_print.h"
#include "fjapp_petsc.h"

using namespace fjlib;


/*!
// Testing ghost cell implementation on a single case
// spherical bubble diffusion, 1 on interface, 0 on the bc
// analytical solution is 0.25(5/r-1) for size=5 case
!*/

int main(int argc, char *args[])
{
	TFJPetsc pet(&argc,&args);
	if (pet.rank==0)
		cout << "using " << pet.size << " computer(s)" << endl;
	bool disp=pet.GetBool("-v",true);
	
	try {
		TFJDiffGC_Coupled solver;
		solver._nlen=pet.GetInt("-size",5);
		solver._npl=pet.GetInt("-npl",2);
		solver._bradius=pet.GetInt("-r",1);
		double tol=pet.GetDouble("-tol",1e-3);
		solver.set_tolerance(tol);
		solver.initialize();
		solver.solve();
		if (disp)
		print_mat("solution",solver.get_cmat(),vmpfText | vmpfColMajor | vmpfColReverse);
		print_var("iters",solver.get_iters());
		save_mat("x.txt",solver.get_cmat());
	}
	catch (const char* s)
	{
		cout << "err: " << s << endl;
	}
	
	if (pet.rank==0)
		cout << "Done." << endl;
	return 0;
}
