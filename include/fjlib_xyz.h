#ifndef FJLIB_XYZ_H
#define FJLIB_XYZ_H

#include "fjlib_cio.h"
#include "fjlib_float.h"

namespace fjlib {

// coordinates direction enumeration type
enum xyz_direction {
	cdX=0,cdY,cdZ
};

// created 3.7.2006
template <class T>
struct xy {
	T x,y;
	xy() {}
	xy(const T& _x, const T& _y) {
		x=_x; y=_y;
	}
	xy(const xy& _xy) {
		if (this!=&_xy) {
			x=_xy.x; y=_xy.y;
		}
	}
};

template <class T>
std::ostream& operator<<(std::ostream& out, const xy<T>& v)
{
	out << v.x << ' ' << v.y;
	return out;
}

template <class T>
std::istream& operator>>(std::istream& in, const xy<T>& v)
{
	in >> v.x >> v.y;
	return in;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

template <class T>
struct xyz {
	T x,y,z;
	xyz() {}
	xyz(const T& _x, const T& _y, const T& _z) {
		x=_x; y=_y; z=_z;
	}
	xyz(const xyz& _xyz) {
		if (this!=&_xyz) {
			x=_xyz.x; y=_xyz.y; z=_xyz.z;
		}
	}
};

}	// end of namespace

#endif
