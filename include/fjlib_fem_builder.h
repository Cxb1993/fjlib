#ifndef FJFEMBuilderH
#define FJFEMBuilderH

// 11.16 allows to build global matrix m and vector b
//		 for mx=b
// 1.3	 size is calculated from coordinates 
//		 instead of stored now
// 1.10	 allows to assign external matrix and vector
//		 through call back function
// 1.14  allows to build nonlinear Jacob and Fx function
//		 storage requires the solution set avariable
//		 and empty 

#include "fjlib_fem_storage.h"
#include "fjlib_function.h"

namespace fjlib {

enum TFJFEMMatVecOptType {
	motInitialize=0,
	motInsert,motAdd,
	motInsertDone,motAddDone,
	motFinalize
};

struct TFJFEMMatVecData {
	int				dim;
	int				i,j;		// -1 means invalid
	float_t			value;
	TFJFEMMatVecOptType
					op;
	void			v(int ii, int jj, float_t vv, 
						TFJFEMMatVecOptType pp)
	{ i=ii; j=jj; value=vv; op=pp; }

	void*			matrix;
	void*			vector;
};

// The purpose of TFJFEMBuilder is to assemble
// global mx=b based on local mx=b of governing eqn,
// integrating BCs if required. Raw data is fed by 
// TFJFEMStorage
class TFJFEMBuilder_Base {
protected:
	matrix_f		mglobal;			// global mx=b
	vector_f		bglobal;
	TFJFEMStorage*
					storage;			// data storage
	// fill element local mx=b into global
	void			merge_mb(const vector_n& ids, 
							const matrix_f& mat,
							const vector_f& vec);
	// fix one node value
	void			diagonal_mb(int id, const float_t v);
public:
	virtual inline
	void			set_storage(TFJFEMStorage* sto)
	{ storage=sto; }
	virtual
	void			initialize()=0;
	virtual	
	void			build()=0;
	matrix_f&		matrix() { return mglobal; }
	vector_f&		vector() { return bglobal; }
protected:
	// call back function to assign matrix
	TFJFEMMatVecData* 
					external_data;
	TFJProcedure	assign_func;
public:
	virtual
	void			set_external(TFJProcedure func,
								TFJFEMMatVecData* data)
	{ assign_func=func; external_data=data; }
	virtual
	void			set_external_dim(int n);
protected:
	bool			is_external_assign()
	{ return !assign_func.empty(); }
	virtual
	void			assign_ex(int i, int j, float_t v, 
						TFJFEMMatVecOptType op)
	{
		external_data->v(i,j,v,op);
		assign_func(external_data);
	}
};

enum TFJFEMBuilderType {
	fbtDefault=0,
	fbtNewton
};

class TFJFEMBuilder: public TFJFEMBuilder_Base {
protected:
	TFJFEMElement
					fem;				// local mx=b
	vector_f		sz;					// size info
	matrix_f		pos;				// nodes pos info
protected:
	TFJFEMEquation	*eqn,				// pde 
					*bc;				// boundary condtions
	TFJFunction		bc_fixed;			// node based bc
private:
	TFJFEMElementData
					data;				// data for pde & bc
public:
	void			set(TFJFEMStorage* sto,
						TFJFEMEquation *_eqn,
						TFJFEMEquation *_bc,
						TFJFunction _fixed)
	{ set_storage(sto); eqn=_eqn; bc=_bc; bc_fixed=_fixed; }
	void			initialize();
	void			build();
public:
	TFJFEMBuilderType
					type;
	virtual inline
	void			set_type(TFJFEMBuilderType t)
	{ type=t; }
	TFJFEMBuilder(): type(fbtDefault) {}
};

}	// end of namespace

#endif