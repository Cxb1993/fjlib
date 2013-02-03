#ifndef FJLIB_TEST_H
#define FJLIB_TEST_H

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
using boost::unit_test::test_suite;

void test_dummy()
{
    BOOST_CHECK(1 == 1); 
}

#define BOOST_TEST_MAIN_MACRO test_suite* \
		init_unit_test_suite(int argc, char* argv[]) { \
			return fjlib_init_test(); \
		} 

#endif
