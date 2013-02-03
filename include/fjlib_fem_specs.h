#ifndef FJFEMElementSPECH
#define FJFEMElementSPECH

// 1.14 add build-in specs for element

#include "fjlib_fem_element.h"

namespace fjlib {

const TFJFEMElementSpec FEMSpec_Line={
					fetLine,
					1,1,2,{-1,1},
					&FEM_SFRect_Linear,
					&FEM_JF,
					&FEM_BFRect };

const TFJFEMElementSpec FEMSpec_Rect={
					fetRect,
					2,1,4,{-1,1},
					&FEM_SFRect_Linear,
					&FEM_JF,
					&FEM_BFRect };

const TFJFEMElementSpec FEMSpec_Triag={
					fetTriangle,
					2,1,3,{0,1},
					&FEM_SFTriag_Linear,
					&FEM_JF,
					&FEM_BFTriag };

}	// end of namespace


#endif