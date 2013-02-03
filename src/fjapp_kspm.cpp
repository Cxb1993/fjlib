#include "fjapp_kspm.h"
#include "fjlib.h"
#include "fjlib_cio.h"

namespace fjlib {

void TFJKSP_Mask::_gen_index_map()
{
	_imat.resize(_mmat->size1(),_mmat->size2());
	reset_matrix(_imat,-1);
	size_t k=0;
	for (size_t i=0; i<_mmat->size1(); i++)
		for (size_t j=0; j<_mmat->size2(); j++)
		{
			if ((*_mmat)(i,j)!=0)	// not known
			{  _imat(i,j)=k; k++; }
		}

	// set total
	_pn=k;
}

void TFJKSP_Mask::initialize()
{
	_gen_index_map();
	TFJKSP_Base::initialize();
}

int TFJKSP_Mask::_uindex(int  m, int n)
{
	if (!between(m,0,(int)_mmat->size1()-1)) return -1;
	if (!between(n,0,(int)_mmat->size2()-1)) return -1;
	return _imat(m,n);
}

void TFJKSP_Mask::_add(PetscInt m, PetscInt n, PetscScalar  v)
{
	if ((v!=0.0) && (n>=0)) {
		_setmat(m,n,v);
//		cout << "add (" << m << "," << n << ") " << v << endl;
	}
}

int	TFJKSP_Mask::build_mat()
{
	PetscInt im,in;
	MatGetOwnershipRange(_pm,&im,&in);
//	cout << im << "-" << in << endl;
	PetscScalar a;
	int	np,nn;
	for (size_t i=0; i<_mmat->size1(); i++)	
		for (size_t j=0; j<_mmat->size2(); j++)
		{
			np=_uindex(i,j);	a=(*_ap)(i,j);
			if ((np<im) || (np>=in)) continue;
			if (a==0.0) throw "ap can't be zero in fjapp_kspm.cpp";
//			_setmat(np,np,a);					// center
			_add(np,np,a);					// center
			_add(np,_uindex(i-1,j),(*_aw)(i,j));// west 
			_add(np,_uindex(i+1,j),(*_ae)(i,j));// east
			_add(np,_uindex(i,j-1),(*_as)(i,j));// south
			_add(np,_uindex(i,j+1),(*_an)(i,j));// north
		}
	// assemble
	_ierr=MatAssemblyBegin(_pm,MAT_FINAL_ASSEMBLY);CHKERRQ2;
	_ierr=MatAssemblyEnd(_pm,MAT_FINAL_ASSEMBLY);CHKERRQ2;

	if (_first_time) {
		_first_time=false;
		MatSetOption(_pm,MAT_NO_NEW_NONZERO_LOCATIONS);
	}

//	MatView(_pm,PETSC_VIEWER_STDOUT_WORLD);

	return 0;
}

int	TFJKSP_Mask::build_vec()
{
	PetscInt im,in;
	VecGetOwnershipRange(_pb,&im,&in);
	PetscScalar a;
	int np;
	for (size_t i=0; i<_mmat->size1(); i++)
		for (size_t j=0; j<_mmat->size2(); j++)
		{
			np=_uindex(i,j);
			if ((np<im) || (np>=in)) continue;
			a=-(*_ab)(i,j);	
			_setvec(np,a);
		}
	// assemble
	_ierr=VecAssemblyBegin(_pb);CHKERRQ2;
	_ierr=VecAssemblyEnd(_pb);CHKERRQ2;
	
//	VecView(_pb,PETSC_VIEWER_STDOUT_WORLD);
	return 0; 
}

int	TFJKSP_Mask::build_x()
{
	Vec x;
	VecScatter sc;
	_ierr=VecScatterCreateToAll(_px,&sc,&x);CHKERRQ2;
	VecScatterBegin(_px,x,INSERT_VALUES,SCATTER_FORWARD,sc);
	VecScatterEnd(_px,x,INSERT_VALUES,SCATTER_FORWARD,sc);

 	PetscScalar *a;
	_ierr=VecGetArray(x,&a);CHKERRQ2;
	int np;
	for (size_t i=0; i<_mmat->size1(); i++)
		for (size_t j=0; j<_mmat->size2(); j++)
		{
			np=_uindex(i,j);
			if (np>=0) (*_ax)(i,j)=a[np];
		}
	_ierr=VecRestoreArray(x,&a);CHKERRQ2;

	VecScatterDestroy(sc); 
	VecDestroy(x);
}

}	// end of namespace
