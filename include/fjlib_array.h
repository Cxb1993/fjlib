#ifndef FJLibArrayH
#define FJLibArrayH

/*!
//	STL vector doesn't define >> and << operator
*/

#include <iostream>
#include <vector>

namespace std {

template <class T>
std::ostream& operator<<(std::ostream& out, std::vector<T>& c)
{
	out << c.size() << endl;
	for (size_t i=0; i<c.size(); i++)
		out << c[i] << endl;
	return out;
}

template <class T>
std::istream& operator>>(std::istream& in, std::vector<T>& c)
{
	int n;
	in >> n;
	if (n<0) return in;
	c.resize(n);
	while (!in.bad() && (n>0))
	{
		in >> c[n-1];
		n--;
	}
	return in;
}

} // end of namespace

#endif
