#ifndef FJLIB_IO
#define FJLIB_IO

#include <iostream>
#include "fjlib_string.h"

using std::cout;
using std::cin;
using std::endl;

namespace fjlib {

template <class T>
void print_var(str_t str, const T& v)
{
	cout << str << ": " << v << endl;
}

template <class T>
void print(const T& v, char space='\n')
{
	cout << v << space;
}

template <class T>
void print_list(str_t str, const T& v, char space='\n')
{
	cout << str << ": " << space;
//	for_each(v.begin(),v.end(),print<T::value_type>);
	for (typename T::const_iterator it=v.begin(); it<v.end(); it++)
		cout << *it << space;
}

}	// end of namespace 

#endif
