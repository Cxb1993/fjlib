#include "fjapp_kspm.h"
#include "fjlib.h"
#include "fjlib_vecmat_print.h"

using namespace fjlib;

int main(int argc, char *args[])
{
	const char help[]="ksp mask test \n\n";
	PetscInitialize(&argc,&args,(char*)0,help);
	PetscMPIInt size,rank;
	CHKERRQ(MPI_Comm_size(PETSC_COMM_WORLD,&size));
	CHKERRQ(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
	cout << "using " << size << " computer(s)" << endl;

	size_t m=8,n=8;
	matrix_n au(m,n);
	matrix_f ap(m,n);
	matrix_f aw(m,n);
	matrix_f ae(m,n);
	matrix_f an(m,n);
	matrix_f as(m,n);
	matrix_f ab(m,n);
	
	// bulk coefficient
	reset_matrix(au,1);		// all unknowns
	reset_matrix(ap,-4.0);
	reset_matrix(aw,1.0);
	reset_matrix(ae,1.0);
	reset_matrix(an,1.0);
	reset_matrix(as,1.0);
	reset_matrix(ab,0.0);

	// diritch b.c
	// left
	for (size_t j=0; j<n; j++)
	{
		as(0,j)=0.0;	an(0,j)=0.0;	aw(0,j)=0.0;	ae(0,j)=0.0;	
		ab(0,j)=4.0*1.0;
		ae(m-1,j)=0.0; 	aw(m-1,j)=0.0;	as(m-1,j)=0.0;	an(m-1,j)=0.0;
	}
	// top and bottom
	for (size_t i=1; i<m-1; i++)
	{
		as(i,0)=0.0;	an(i,0)=4.0;	aw(i,0)=0.0;	ae(i,0)=0.0;
		an(i,n-1)=0.0;	as(i,n-1)=4.0;	aw(i,n-1)=0.0;	ae(i,n-1)=0.0;
	}
/*	
	print_mat("ap",ap,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("aw",aw,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("ae",ae,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("as",as,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("an",an,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("ab",ab,vmpfText | vmpfColMajor | vmpfColReverse);
*/
	// assemebly and solve
	matrix_f ax(m,n);		// solution
	TFJKSP_Mask solver;
	solver.set_data(&ap,&aw,&ae,&as,&an,&ab,&ax);
	solver.set_mask(&au);
	solver.set_tolerance(1e-12);
	solver.initialize();
//	print_mat("imat",solver.get_imat(),vmpfText | vmpfColMajor | vmpfColReverse);
	solver.solve();
	print_mat("ax",ax,vmpfText | vmpfColMajor | vmpfColReverse);

	solver.finalize();
	PetscFinalize();
	cout << "Done." << endl;
	return 0;
}
