#ifndef FJLIBVecMatPrintH
#define FJLIBVecMatPrintH

#include "fjlib_vecmat.h"
//#include <iostream>
#include "fjlib_cio.h"
#include "fjlib_string.h"
#include <boost/format.hpp>
#include "fjlib_io.h"
/*!
//	1.5 replace char space by spacing using boost format
!*/

namespace fjlib {

enum TFJVecMatPrintFormat {
	vmpfText=0x000001L,				///<- default if not set
	vmpfNoHeader=0x000002L,			///<- header info, ex dimension
	vmpfReverse=0x000004L,			///<- reverse the order for matrix
	vmpfColMajor=0x000008L,			///<- major, default row major
	vmpfRowReverse=0x000010L,		///<- reverse the order of rows
	vmpfColReverse=0x000020L,		///<- reverse the order for cols
	vmpfMultiplyVec=0x000040L,		///<- read matrix vector by vector
	vmpfSingleVec=0x000080L			///<- read single vector into matrix
};

template <class T>
class TFJVecPrint {
public:
	typedef boost::numeric::ublas::vector<T> vector_type;
protected:
	vector_type*	vec;
	long			options;
	char			space;
public:
	TFJVecPrint(): space(' ') {}
    TFJVecPrint(vector_type* v, long _options=vmpfText):
		  vec(v), options(_options), space(' ') {}
	void			set(vector_type* v, long _options=vmpfText)
	{ vec=v; options=_options; }
	void			set_vector(vector_type* v)
	{ vec=v; }
	void			set_options(long _options) 
	{ options=_options; }
	void			set_space(char _space)
	{ space=_space; }

	template <class T2> friend
	std::ostream& operator<<(std::ostream& os, const TFJVecPrint<T2>& v);
	template <class T2> friend
	std::istream& operator>>(std::istream& in, TFJVecPrint<T2>& v);
};

template <class T2>
std::ostream& operator<<(std::ostream& os, const TFJVecPrint<T2>& v)
{
	const typename TFJVecPrint<T2>::vector_type& vec=*(v.vec);
	if (!(vmpfText & v.options)) { os << vec; return os; }

	size_t n=vec.size();
	if (!(vmpfNoHeader & v.options))
		os << n << std::endl;

	if (n<1) return os;
	if (!(vmpfReverse & v.options))
		for (size_t i=0; i<n; i++)
		{
			os << vec[i];
			if (vmpfColMajor & v.options) os << std::endl;
			else os << v.space;
		}
	else
		for (int i=n-1; i>=0; i--)
		{
			os << vec[i];
			if (vmpfColMajor & v.options) os << std::endl;
			else os << v.space;
		}
	return os;
}

template <class T2>
std::istream& operator>>(std::istream& in, TFJVecPrint<T2>& v)
{
	typename TFJVecPrint<T2>::vector_type& vec=*(v.vec);
	size_t s;
	in>>s;
	vec.resize(s);
	size_t i=0;
	while ((i<s) && (in))
	{
		in>>vec[i];
	    ++i;
	}
	return in;
}

typedef TFJVecPrint<size_t> vprint_sz;
typedef TFJVecPrint<int> vprint_n;
typedef TFJVecPrint<float_t> vprint_f;
typedef TFJVecPrint<bool> vprint_b;

template <class T>
void print_vec(boost::numeric::ublas::vector<T>& v, long options=vmpfText)
{
	TFJVecPrint<T> p(&v,options);
	cout << p << endl;
}

template <class T>
void print_vec(str_t title, boost::numeric::ublas::vector<T>& v, long options=vmpfText)
{
	cout << title << ": " << endl;
	print_vec(v,options);
}

template <class T>
void save_vec(str_t fn, boost::numeric::ublas::vector<T>& v, long options=vmpfText)
{
	TFJVecPrint<T> p(&v,options);
	quicksave(fn.c_str(),p);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class T>
class TFJMatPrint { 
public:
	typedef boost::numeric::ublas::matrix<T> matrix_type;
protected:
	matrix_type*	mat;
	long			options;
//	char			space;
	size_t			space;
public:
	TFJMatPrint(): options(vmpfText), space(8) {}
    TFJMatPrint(matrix_type* v, long _options=vmpfText):
		  mat(v), options(_options), space(8) {}
	void			set(matrix_type* v, long _options=vmpfText)
	{ mat=v; options=_options; }
	void			set_matrix(matrix_type* v)
	{ mat=v; }
	void			set_options(long _options) 
	{ options=_options; }
	void			set_spacing(size_t _space)
	{ space=_space; }
	long			get_options()
	{ return options; }

	template <class T2> friend
	std::ostream& operator<<(std::ostream& os, const TFJMatPrint<T2>& v);
	template <class T2> friend
	std::istream& operator>>(std::istream& in, TFJMatPrint<T2>& v);
};

template <class T2>
std::ostream& operator<<(std::ostream& os, const TFJMatPrint<T2>& v)
{
	const typename TFJMatPrint<T2>::matrix_type& mat=*(v.mat);
	if (!(v.options & vmpfText)) { os << mat; return os; }

	size_t m=mat.size1(), n=mat.size2();
	if (!(vmpfNoHeader & v.options))
		os << m << ' ' << n << std::endl;

//	str_t fm="%|"+ to_string(v.space) + "|";
	str_t fm="% "+ to_string(v.space) + "f";
	size_t ii,jj;
	if (!(vmpfColMajor & v.options))
	{
		for (size_t i=0; i<m; i++)
		{
			ii=i; if (vmpfRowReverse & v.options) ii=m-1-i;
			for (size_t j=0; j<n; j++)
			{
				jj=j; if (vmpfColReverse & v.options) jj=n-1-j;
//		os << mat(ii,jj) << v.space;
				os << boost::format(fm) % mat(ii,jj) << " ";
			}
			os << std::endl;
		}
	}
	else
	{
		for (size_t j=0; j<n; j++)
		{
			jj=j; if (vmpfColReverse & v.options) jj=n-1-j;
			for (size_t i=0; i<m; i++)
			{
				ii=i; if (vmpfRowReverse & v.options) ii=m-1-i;
//				os << mat(ii,jj) << v.space;
				os << boost::format(fm) % mat(ii,jj) << " ";		
			}
			os << std::endl;
		}
	}
	return os;	
}

template <class T2>
std::istream& operator>>(std::istream& in, TFJMatPrint<T2>& v)
{
	typename TFJMatPrint<T2>::matrix_type& mat=*(v.mat);
	if (vmpfSingleVec & v.options)
	{
		typedef boost::numeric::ublas::vector<T2> vector_type;
		vector_type vec;
		TFJVecPrint<T2> vec_ptr(&vec);
		in >> vec_ptr;
		mat.resize(1,vec.size());
		for (size_t i=0; i<vec.size(); i++)
			mat(0,i)=vec[i];
		return in;
	}
	if (vmpfMultiplyVec & v.options)
	{
		size_t m;
		in >> m;
		// load first vector
		typedef boost::numeric::ublas::vector<T2> vector_type;
		vector_type vec;
		TFJVecPrint<T2> vec_ptr(&vec);
		in >> vec_ptr;
		size_t n=vec.size(); 

		mat.resize(m,n);
		int i=0;
		while ((i<m) && in)
		{
			// copy vec
			for (size_t j=0; j<n; j++)
				mat(i,j)=vec[j];
			i++;
			if (i<m)
				in >> vec_ptr;
		}
		return in;
	}

	// default loading 
		size_t m,n;
		in>>m>>n;
		mat.resize(m,n);
		size_t i=0, j=0;
		while ((i<mat.size1()) && (in))
		{
			in>>mat(i,j);
			++j; 
			if (j==mat.size2()) 
			{ j=0; ++i; }
		}

	return in;
}


typedef TFJMatPrint<size_t>		mprint_sz;
typedef TFJMatPrint<int>		mprint_n;
typedef TFJMatPrint<float_t>	mprint_f;
typedef TFJMatPrint<bool>		mprint_b;

template <class T>
void print_mat(boost::numeric::ublas::matrix<T>& m, long options=vmpfText)
{
	TFJMatPrint<T> p(&m,options);
	cout << p << endl;
}

template <class T>
void print_mat(str_t title, boost::numeric::ublas::matrix<T>& m, long options=vmpfText)
{
	cout << title << ": " << endl;
	print_mat(m,options);
}

template <class T>
void save_mat(str_t fn, boost::numeric::ublas::matrix<T>& m, long options=vmpfText)
{
	TFJMatPrint<T> ptr(&m,options);
	quicksave(fn.c_str(),ptr);
}


} // end of namespace

#endif
