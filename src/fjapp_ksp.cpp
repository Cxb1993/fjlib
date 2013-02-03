#include "fjapp_ksp.h"
#include "fjlib_cio.h"

namespace fjlib {

int	TFJKSP_Base::create_ksp()
{
	_ierr=KSPCreate(PETSC_COMM_WORLD,&_ksp);CHKERRQ2;
	_ierr=KSPSetFromOptions(_ksp);CHKERRQ2;
	_first_time=true;
}

int	TFJKSP_Base::destroy_ksp()
{
	_ierr=KSPDestroy(_ksp);CHKERRQ2;
}

int TFJKSP_Base::create_matvec()
{
//	if (_created) destroy_matvec();
	if (_created) return 0;

	// create b and x vector
	_ierr=VecCreate(PETSC_COMM_WORLD,&_pb);CHKERRQ2;
	_ierr=VecSetSizes(_pb,PETSC_DECIDE,_pn);CHKERRQ2;
	_ierr=VecSetFromOptions(_pb);CHKERRQ2;
	_ierr=VecDuplicate(_pb,&_px);CHKERRQ2;

	// create m matrix, changes made for 2.3.0
	_ierr=MatCreate(PETSC_COMM_WORLD,&_pm); CHKERRQ2;
	_ierr=MatSetSizes(_pm,PETSC_DECIDE,PETSC_DECIDE,_pn,_pn);CHKERRQ2;
	_ierr=MatSetFromOptions(_pm);CHKERRQ2;
	
	_created=true;
}

int TFJKSP_Base::destroy_matvec()
{
	if (!_created) return 0;
	_ierr=VecDestroy(_pb);CHKERRQ2;
	_ierr=VecDestroy(_px);CHKERRQ2;
	_ierr=MatDestroy(_pm);CHKERRQ2;
	_created=false;
	_first_time=true;
}

int	TFJKSP_Base::solve_matvec()
{
	// set options and preconditioner
	_ierr=KSPSetOperators(_ksp,_pm,_pm,
				SAME_NONZERO_PATTERN);CHKERRQ2;
//	_ierr=KSPSetOperators(_ksp,_pm,_pm,
//				DIFFERENT_NONZERO_PATTERN);CHKERRQ2;
	_ierr=KSPSetTolerances(_ksp,_tol,PETSC_DEFAULT,
				PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ2;
	// solve it
	_ierr=KSPSolve(_ksp,_pb,_px);CHKERRQ2;
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

	_ierr=KSPGetIterationNumber(_ksp,&_iters);CHKERRQ2;
}

/*
int	TFJKSP_Base::clear_matvec()
{
	PetscScalar a=0.0;
	VecSet(pb,a);
//	_ierr=VecSet(&a,pb);CHKERRQ2;	// changs made for 2.3.0

//	_ierr=MatZeroEntries(pm);CHKERRQ2;
}
*/

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

int TFJKSP::build_vec()
{
	PetscScalar a;
	int	np;
	for (size_t i=0; i<_nx; i++)
		for (size_t j=0; j<_ny; j++)
		{
			np=j*_nx+i;
			// Update vector
			a=-(*_ab)(i+1,j+1);
			_setvec(np,a);CHKERRQ2;
		}

	// assemble
	_ierr=VecAssemblyBegin(_pb);CHKERRQ2;
	_ierr=VecAssemblyEnd(_pb);CHKERRQ2;
}

int TFJKSP::build_mat()
{
	PetscInt im,in;
	MatGetOwnershipRange(_pm,&im,&in);
	PetscScalar a;
	int	np;
	for (size_t i=0; i<_nx; i++)	
		for (size_t j=0; j<_ny; j++)
		{
			np=j*_nx+i;
			if ((np<im) or (np>=in)) continue;

			// Update matrix
				// center
			a=(*_ap)(i+1,j+1);
			_setmat(np,np,a);CHKERRQ2;
				// west
			if (i>0) { 
				a=(*_aw)(i+1,j+1);
				_setmat(np,np-1,a);CHKERRQ2;
			}
				// east
			if (i<_nx-1) {
				a=(*_ae)(i+1,j+1);
				_setmat(np,np+1,a);CHKERRQ2;
			}
				// south
			if (j>0) {
				a=(*_as)(i+1,j+1);
				_setmat(np,np-_nx,a);CHKERRQ2;
			}
				// north
			if (j<_ny-1) {
				a=(*_an)(i+1,j+1);
				_setmat(np,np+_nx,a);CHKERRQ2;
			}
		}
	// assemble
	_ierr=MatAssemblyBegin(_pm,MAT_FINAL_ASSEMBLY);CHKERRQ2;
	_ierr=MatAssemblyEnd(_pm,MAT_FINAL_ASSEMBLY);CHKERRQ2;

	if (_first_time) {
		_first_time=false;
		MatSetOption(_pm,MAT_NO_NEW_NONZERO_LOCATIONS);
	}
}


int TFJKSP::build_x()
{
	// setup a local vec to scatter the data
	Vec x;
	VecScatter sc;
	_ierr=VecScatterCreateToAll(_px,&sc,&x);CHKERRQ2;
	VecScatterBegin(_px,x,INSERT_VALUES,SCATTER_FORWARD,sc);
	VecScatterEnd(_px,x,INSERT_VALUES,SCATTER_FORWARD,sc);

 	PetscScalar *a;
	_ierr=VecGetArray(x,&a);CHKERRQ2;
	for (size_t i=0; i<_nx; i++)
		for (size_t j=0; j<_ny; j++)
			(*_ax)(i+1,j+1)=a[j*_nx+i];
	_ierr=VecRestoreArray(x,&a);CHKERRQ2;

	VecScatterDestroy(sc); 
	VecDestroy(x);
}

}
