#ifndef FJLIB_BUBBLE_NEEDLE_NSH
#define FJLIB_BUBBLE_NEEDLE_NSH

#include "fjapp_bubbleNS.h"

namespace fjlib {

/// simple wall position type
enum TFJWallPosType {
	bdtS=0,bdtE,bdtN,bdtW
};

/// boundary condition type
enum TFJBCType {
	bctVoid,			// None, excluded
	bctDirichet,		// u=0
	bctNeumann,			// du=0
	bctUnknown			// user supplied
};

/// bc maker info data
struct TFJBCMarkerInfoType {
	int i,j,k;	
	TFJWallPosType pos;
	TFJBCType type;
	float_t value[3];
};

/*!
// extension work of bubbleNS
// needle or sold walls are added with non-slip bc, therefore
// bcMarker are introduced to correct the equation formulation
// near/on wall regions
//
// Cation:
// TFJBCMarkerInfoType pos uses internal coordinates axis, sptO
//
// 7.20	created
!*/
class TFJBubbleNeedleNS_Solver: public TFJBubbleNS_Solver {
private:
	/// used for mark the boundary and type
	vector_f	m_BCMarkerList;
	void		_set_dirichet_value(matrix_f& ap, 
					matrix_f& aw, matrix_f& ae, 
					matrix_f& as, matrix_f& an, 
					matrix_f& su, 
					size_t i, size_t j, float_t v);
	void		_set_neumann(matrix_f& ap, 
					matrix_f& aw, matrix_f& ae, 
					matrix_f& as, matrix_f& an, 
					matrix_f& su, 
					size_t i, size_t j, TFJWallPosType wp);
protected:
	void		prepare_bc(staggered_type st);
	void		set_dirichet(staggered_type st,
						size_t i, size_t j, float_t v);
public:
	void		clear_BCMarker() { m_BCMarkerList.resize(0); }
};

}	// end of namespace

#endif
