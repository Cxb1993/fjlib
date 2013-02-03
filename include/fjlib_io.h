#ifndef FJLIBIOH
#define FJLIBIOH

#include <iostream>
#include <fstream>

namespace fjlib {

template <class T>
void quicksave(const char* fn, const T& v, bool bf=false)
{
	std::ofstream of;
	if (!bf)
		of.open(fn,std::ios::out);
	else
		of.open(fn,std::ios::out | std::ios::binary);

	if (of.good())
		of << v;
	of.close();
}

template <class T>
void quickload(const char* fn, T& v, bool bf=false)
{
	std::ifstream in;
	if (!bf)
		in.open(fn,std::ios::in);
	else
		in.open(fn,std::ios::in | std::ios::binary);

	if (in.good())
		in >> v;
	in.close();
}

/*
template <class T>
void quicksave2(const char* fn, T& v)
{
	std::ofstream of;
	of.open(fn,std::ios::out);
	if (of.good())
		v.ToFile(of);
	of.close();
}

template <class T>
void quickload2(const char* fn, T& v)
{
	std::ifstream ifs;
	ifs.open(fn,std::ios::in);
	if (ifs.good())
		v.LoadFile(ifs);
	ifs.close();
}
*/

}  // end of namespace


//---------------------------------------------------------------------------
#endif
