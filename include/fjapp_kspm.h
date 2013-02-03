#ifndef FJAPP_KSP_MASKED_H
#define FJAPP_KSP_MASKED_H

// Parallel solver for 5-nodes 2D linear system
// Mask is supported by using a mask map matrix

// 1.10.06 created. supposed to be derived from TFJKSP
//			by adding a node mapping routine _uindex(m,n)
//			but can't do it before fixing matrix range 
//			for now it's just a seperate class from TFJKSP_Base
// 2.27.06 fix a bug, _set_vec(a)
// Note: for the mask, if the value equals to non-zero, the cell is treated as unknown
#include "fjapp_ksp.h"

namespace fjlib {

class TFJKSP_Mask: public TFJKSP_Base {
private:
	void 		_add(PetscInt m, PetscInt n, PetscScalar  v);
protected:
	matrix_n	*_mmat;		// unknown map matrixp
	matrix_n	_imat;		// index matrix

	// create one-one map between unknown map and solver
	virtual
	void		_gen_index_map();
	virtual
	int			_uindex(int m, int n);
public:
	TFJKSP_Mask(): TFJKSP_Base() {}
	void		set_mask(matrix_n* mat)
	{ _mmat=mat; }
	void 		initialize();
	int			build_vec();
	int			build_mat();
	int			build_x();
	matrix_n&	get_imat() { return _imat; }
};

}	// end of namespace



#endif
