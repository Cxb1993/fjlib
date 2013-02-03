#ifndef FJLIBKSPH
#define FJLIBKSPH

// Interface between fjlib and petsc ksp solver
// for structured uniform mesh only

// 2D 5 nodes paraell solver for linear systems

/*!
// 1.10.06 revised to include TFJKSP_Base class for further derivation of TFJKSP_Masked
//			the class is one time job, solve different pattern matrix needs to create another of this class
// 1.26.06 changed the finalize mechanism back to c++ compatible, 
//			wrap the code in between {} before petsc_finalize()
!*/

#include "fjlib_vecmat.h"
#include "petscksp.h"

#ifndef _Petsc_KSP
#define _Petsc_KSP
#endif

namespace fjlib {

#define CHKERRQ2 CHKERRQ(_ierr);

/*!
//	Caution: call finalize() before end mpi session
//			for now I don't know how to make it fully c++
//	2.27.06 added setup(), call it for solving another matrix
//			with different structure, ksp bug!!
!*/
class TFJKSP_Base {
protected:
	PetscErrorCode 
				_ierr;			// error storage, reserved
	bool		_created;		// check if matvec created
//	bool		_ksp_created;	// check if ksp created
	bool		_first_time;	// check if first time
	PetscInt	_pn;			// unknown count;
	Vec			_pb,_px;		// rhs and unknowns vectors
	Mat			_pm;			// linear matrix
	PetscScalar	_tol;			// tolerance
	KSP			_ksp;			// solver object
	PCType		_pc;			// preconditioner object
	int			_iters;			// final iteration count

	/// 2D 5-nodes cooefficient layout
	matrix_f	*_ax,*_ap,*_ae,*_aw,*_an,*_as,*_ab;
protected:
	/// Set rhs vector value to v at index m
	inline
	void		_setvec(int m, const PetscScalar& v)
	{ _ierr=VecSetValue(_pb,m,v,INSERT_VALUES); }
	/// Set global matrix value to v at index m,n
	inline
	void		_setmat(int m, int n, const PetscScalar& v)
	{ _ierr=MatSetValue(_pm,m,n,v,INSERT_VALUES); }
public:
	TFJKSP_Base(): _created(false), _tol(1e-5)
	{ create_ksp(); }
	~TFJKSP_Base() { finalize(); }
	/// have to call this before MPI Finalize() call
	/// nothing is allowed after this call
	void		set_data(matrix_f *qp,
						matrix_f *qw, matrix_f *qe,
						matrix_f *qs, matrix_f *qn,
						matrix_f *qb, matrix_f *qx)
	{
		_ap=qp; _aw=qw; _ae=qe; _as=qs;
		_an=qn; _ab=qb; _ax=qx;
	}
	void		set_tolerance(float_t t) { _tol=t; }
protected:
	void		finalize() 
	{ destroy_matvec(); destroy_ksp(); }
	int			create_ksp();
	int			destroy_ksp();
	int			create_matvec();
//	int			clear_matvec();
	/// Build rhs vector
	virtual 
	int			build_vec()=0;
	/// Build global matrix
	virtual
	int			build_mat()=0;
	int			destroy_matvec();
	int			solve_matvec();
	/// Import result into 2D format
	virtual
	int			build_x()=0;
	/// Real unknonw index in the global matrix
public:
	/// Return total unknown count
	size_t		count() { return _pn; }
	/// Initialize the solver, including setup global matrix
	virtual
	void		initialize() 
	{ create_matvec(); build_mat(); }
	/// Solve the final matrix and update solution
	virtual
	void		solve()
	{
		build_vec(); solve_matvec(); build_x();
	}
	/// Returns the total number of interation spent for the last solving 
	const int	niters() const { return _iters; }

	void 		reset()
	{ finalize(); create_ksp(); }
};	// end of TFJKSP_Base


/*!
//
// 1.11.06 revised to be derived from TFJKSP_Base
//			need further improvement to include matrix range
//			right now, the outter edge is excluded.
!*/
class TFJKSP: public TFJKSP_Base {
protected:
	size_t		_nx,_ny;
public:
	TFJKSP(): TFJKSP_Base() {}
	TFJKSP(size_t m, size_t n): TFJKSP_Base()
	{ set_size(m-2,n-2); } 
	void		set_size(size_t m, size_t n)
	{ _nx=m-2; _ny=n-2; _pn=_nx*_ny; 	}	
protected:
	int			build_vec();
	int			build_mat();
	int			build_x();
};	// end of TFJKSP

}	// end of namespace

#endif


