#ifndef FJAPP_PETSC
#define FJAPP_PETSC

#include "petsc.h"

namespace fjlib {

/*!
// created 1.26.06 to faciliate creating petsc app
// wrap other petsc objects within {} or catch {} to ensure
// they are automatically managed
!*/

class TFJPetsc {
public:
	PetscMPIInt size,rank;
protected:
	int	_petsc_init(int* argc, char **args[])
	{
		const char help[]="petsc";
		PetscInitialize(argc,args,(char*)0,help);
		CHKERRQ(MPI_Comm_size(PETSC_COMM_WORLD,&size));
		CHKERRQ(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
		_binit=true;
	}
private:
	bool 	_binit;
public:
	TFJPetsc():_binit(false) {}
	TFJPetsc(int* argc, char** args[]):_binit(false)
	{ Init(argc,args); }
	~TFJPetsc() 
	{ if (_binit) PetscFinalize(); }
	bool 	Init(int *argc, char** args[]) 
	{ if (!_binit) _petsc_init(argc,args); }
	bool	GetBool(char *str,bool v)
	{
		PetscTruth t=(PetscTruth)v;
		PetscOptionsGetTruth(PETSC_NULL,str,&t,PETSC_NULL);
		return t;
	}
	int		GetInt(char *str,int v)
	{
		PetscInt t=v;
		PetscOptionsGetInt(PETSC_NULL,str,&t,PETSC_NULL);
		return t;
	}
	double	GetDouble(char *str,double v)
	{
		PetscScalar t=v;
		PetscOptionsGetScalar(PETSC_NULL,str,&t,PETSC_NULL);
		return t;
	}
};	 // end of class


}	// end of namespace

#endif
