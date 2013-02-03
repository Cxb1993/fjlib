#ifndef FJLIB_POINT_H
#define FJLIB_POINT_H

#include "fjlib_xyz.h"

namespace fjlib {

typedef xy<size_t> point2_sz;
typedef xy<int> point2_n;
typedef xy<float_t> point2_f;

enum CoordEnumXYZ { ceX=0,ceY,ceZ };

}	// end of namespace
#endif
