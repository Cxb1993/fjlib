#ifndef FJFEMElementH
#define FJFEMElementH

// 11.15 allows to build local matrix m and vector b
//		 for mx=b
// 11.16 element size is stand-alone parameter
// 1.3	 merge element and bc data object
//		 difference between them is that bc uses integration
//		 with fixed coordinates TFJIntegration_NDim
// 1.6	 add element nodes postion, more general approach
//		 set_element() allows access to TFJFEMElementStore
//		 element size option still supported
// 1.14	 add non-linear solver capability
//		 set_values() supplies extra storage for nodes value
//		 set eqn.m=Jacob, eqn.b=Fx for newton method

#include "fjlib_float.h"
#include "fjlib_vector.h"
#include "fjlib_matrix.h"
#include "fjlib_function.h"
#include "fjlib_integration_ndim.h"

namespace fjlib {

enum TFJFEMElementName {
	fetLine, fetRect, fetTriangle
};

enum TFJFEMElementOrder {
	feoLinear=1, feoQuadratic, feoCubic
};

struct TFJFEMElementSpec {
	TFJFEMElementName
				name;		
	int			dim;		// dimension
	int			order;		// order
	int			nnodes;		// number of nodes

	float_t		range[2];	// local coordinates range
	TFJProcedure
				shape;		// shape function
	TFJProcedure
				jacob;		// size transformation
	TFJProcedure 
				binfo;		// boundary integration info
};

// type for both equation and boundary
struct TFJFEMEquation {
	TFJFunction1p
				m;			// for stiff matrix m
	TFJFunction1p
				b;			// for b in mx=b
	vector_f	params;		// parameters of the eqn
};

struct TFJFEMElementData: public TFJIntegrate_Data {
	TFJFEMElementSpec*
				spec;
	TFJFEMEquation*
				eqn;
	matrix_f*	npos;		// nodes position
	vector_f*	size;		// size of the element

	// it's up to the jacob and shape to use
	// either npos or size, or both of them

	matrix_f	jacob;		// hold value from shape function
	matrix_f	shape;		// hold jacob matrix
	int			index[2];	// identification of node-node

	matrix_f*	values;		// values or properties 
							// which used for nonlinear system

	void*		element;	// reserved for pointer of 
							// TFJFEMElementStore if needed
	int			bc;			// reserved for boundary cond							
	matrix_f	bfixed;		// hold value for compact bc info
};

// TFJFEMElement has the functionality of building
// local mx=b based on either governing eqn or BC eqn.
// It doesn't store the element information
class TFJFEMElement {
protected:
	TFJFEMElementData*		// convenience for derived data type
				data;
	TFJFEMElementSpec*
				spec;		// specification
	TFJFEMEquation*
				eqn;
	matrix_f	mlocal;
	vector_f	blocal;
private:
	TFJIntegrate_NDim
				integ;
	TFJIntegrate_GaussQuad 
				gq[3];
protected:
	virtual
	void		initialize_integration();
	virtual
	void		initialize_data();
protected:
	void*		_element;
	int			_bc;
	matrix_f*	_values;
public:
	virtual inline
	void		set_data(TFJFEMElementData* d)
	{ data=d; }
	virtual inline
	void		set(TFJFEMElementSpec* s, 
					TFJFEMEquation* e)
	{ spec=s; eqn=e; }
	virtual inline
	void		set_element(void* el)
	{ _element=el; }
	virtual inline
	void		set_bc(int bc)
	{ _bc=bc; }
	virtual inline
	void		set_values(matrix_f *val)
	{ _values=val; }
	virtual
//	void		build(vector_f* size);
	void		build(matrix_f* pos, vector_f* size);
	matrix_f&	matrix() { return mlocal; }
	vector_f&	vector() { return blocal; }
private:		// for fixed integration
	vector_b	fixed;
	vector_f	fixed_value;
protected:
	virtual
	void		initialize_bc(int bc);
	inline
	bool		is_bc() { return (_bc>=0); }
};

// General Jocab 
void FEM_JF(void* obj);

// 1D linear element
void FEM_JF1D(void* obj);
void FEM_SF1D_Linear(void* obj);

// 2D rectangular linear element
void FEM_JFRect(void* obj);
void FEM_SFRect_Linear(void* obj);
void FEM_BFRect(void* obj);

// 2D triangular linear element
void FEM_SFTriag_Linear(void* obj);
void FEM_BFTriag(void* obj);


}	// end of namespace

#endif