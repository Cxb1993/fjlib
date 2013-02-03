#ifndef FJLibFloatH
#define FJLibFloatH

#include <boost/limits.hpp>

namespace fjlib {

typedef double float_t;

#define FJLIB_FLOAT_EPS std::numeric_limits<float_t>::epsilon()

}

#endif

