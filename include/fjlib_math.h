#ifndef TFJLIB_MATH
#define TFJLIB_MATH

#include <cmath>

namespace fjlib {

float round(const float &number, const int num_digits)
{
	float doComplete5i, doComplete5(number * powf(10.0f, (float) (num_digits + 1)));
	
	if(number < 0.0f)
		doComplete5 -= 5.0f;
	else
		doComplete5 += 5.0f;
	
	doComplete5 /= 10.0f;
	modff(doComplete5, &doComplete5i);
	
	return doComplete5i / powf(10.0f, (float) num_digits);
}


}	// end of namespace

#endif