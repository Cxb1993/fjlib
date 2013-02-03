#ifndef TFJLIBFunctionH
#define TFJLIBFunctionH

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "fjlib_float.h"

namespace fjlib {

typedef boost::function<void(void *) >		TFJProcedure;
typedef boost::function<float_t(void *) >	TFJFunction;
typedef boost::function<float_t(void *,float_t *) >
											TFJFunction1p;
typedef boost::function<float_t(void *, float_t *, float_t *) >
											TFJFunction2p;
typedef boost::function<float_t(void *, vector_f *) >
											TFJFunctionNp;
typedef boost::function<void(void *, float_t *, vector_f *) >
											TFJProcIn1OutNp;
}	// end of namespace

#endif