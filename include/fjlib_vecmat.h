#ifndef FJLIB_VECMAT_H
#define FJLIB_VECMAT_H

// seperate ublas & stl matrix 3/19/2006
// both works, just define USE_STL_MATRIX to use standard

#ifndef USE_STL_MATRIX
#include "fjlib_vecmat_ublas.h"
#else
#include "fjlib_vecmat_stl.h"
#endif

#include "fjlib_float.h"
#include <stdarg.h>
namespace fjlib {

// typedefs
typedef VECTOR<size_t> vector_sz;
typedef VECTOR<int> vector_n;
typedef VECTOR<float_t> vector_f;
typedef VECTOR<bool> vector_b;

typedef MATRIX<size_t> matrix_sz;
typedef MATRIX<int> matrix_n;
typedef MATRIX<bool> matrix_b;
typedef MATRIX<float_t> matrix_f;

//////////////////////////////////////////////////////////
/////////////////// Vector Utilities /////////////////////
//////////////////////////////////////////////////////////

template <class T> inline
void resize(VECTOR<T>& v1, VECTOR<T>& v2)
{
	v1.resize(v2.size());
}

template <class T> inline
void push_back(VECTOR<T>& vec, const T& v)
{
	size_t n=vec.size();
	vec.resize(n+1);
	vec[n]=v;
}

template <class T> inline
void pop_back(VECTOR<T>& vec)
{
	int n=vec.size();
	if (n>0) vec.resize(n-1);
}

template <class T> 
void insert(VECTOR<T>& vec, size_t pos, const T& v)
{
	int n=vec.size()+1;
	if (pos>=n) throw;
	vec.resize(n);
	if (n==1) { vec[0]=v; return; }
	for (int i=n-1; i>pos; i--)
		vec[i]=vec[i-1];
	vec[pos]=v;
}

template <class T>
void append(VECTOR<T>& vec, const VECTOR<T>& v)
{
	if (v.size()<1) return;
	int oldn=vec.size();
	int newn=oldn+v.size();
	vec.resize(newn);
	for (size_t i=0; i<v.size(); i++)
		vec[i+oldn]=v[i];
}

// slice v(p0,p1) into v1
template <class T> inline
void slice(const VECTOR<T>& v,
			VECTOR<T>& v1, size_t p0, size_t p1)
{
	v1.resize(p1-p0+1);
	copy(v.begin()+p0,v.begin()+p1+1,v1.begin());
}

template <class T> 
void split(const VECTOR<T>& v, 
		   VECTOR<T>& v1, VECTOR<T>& v2,
		   size_t v1_count, bool joint_copy=false)
{
	if (v1_count==0)
	{ v1.resize(0); v2=v; return; }

	if (v1_count>v.size()) throw "v1count>size";

	slice(v,v1,0,v1_count-1);
	slice(v,v2,v1_count-joint_copy,v.size()-1);
}

template <class T>
void erase(VECTOR<T>& v, size_t p0, size_t p1)
{
	copy(v.begin()+p1+1,v.end(),v.begin()+p0);
	v.resize(v.size()-(p1-p0+1));
}

template <class T>
void reset_vector(VECTOR<T> &vec,const T& v=0)
{
	if (vec.size()==0) return;
	for (size_t i=0; i<vec.size(); i++)
		vec[i]=v;
}

// don't call this directly, it dones't have safe type check 
// call set_vector instead
template <class T>
void set_vec(VECTOR<T>& vec, int n, T v, ...)
{
	if (n<1) throw "wrong size in set_vec in fjlib_vecmat.h";
	size_t nn=n;
    vec.resize(nn);
	va_list ap;
	va_start(ap,v);
	for (size_t i=0; i<nn; i++)
	{
		vec[i]=v;
		v = va_arg(ap, T);
	}
	va_end(ap);	
}

template <class T>
void set_vector(VECTOR<T> &v, const T& v0)
{v.resize(1); v[0]=v0;}
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1)
{v.resize(2); v[0]=v0; v[1]=v1; }
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2)
{v.resize(3); v[0]=v0; v[1]=v1; v[2]=v2; }
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2, const T& v3)
{set_vec(v,4,v0,v1,v2,v3);}
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2, const T& v3, const T& v4)
{set_vec(v,5,v0,v1,v2,v3,v4);}
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5)
{set_vec(v,6,v0,v1,v2,v3,v4,v5);}
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6)
{set_vec(v,7,v0,v1,v2,v3,v4,v5,v6);}
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7)
{set_vec(v,8,v0,v1,v2,v3,v4,v5,v6,v7);}
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7, const T& v8)
{set_vec(v,9,v0,v1,v2,v3,v4,v5,v6,v7,v8);}
template <class T>
void set_vector(VECTOR<T> &v, const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7, const T& v8, const T& v9)
{set_vec(v,10,v0,v1,v2,v3,v4,v5,v6,v7,v8,v9);}

//////////////////////////////////////////////////////////
/////////////////// Matrix Utilities /////////////////////
//////////////////////////////////////////////////////////

enum TFJRowCol {
	rcRow=0,
	rcCol
};

template <class T> inline
void resize(MATRIX<T>& m1, MATRIX<T>& m2)
{
	m1.resize(m2.size1(),m2.size2());
}

template <class T> 
void reset_matrix(MATRIX<T> &mat,const T& v=0)
{
	if ((mat.size1()==0) || (mat.size2()==0)) return;
	for (size_t i=0; i<mat.size1(); i++)
		for (size_t j=0; j<mat.size2(); j++)
			mat(i,j)=v;
}

// make a m1=m((tlr,tlc)-(brr,brc))
template <class T>
void slice_matrix(const MATRIX<T>& m,
		   MATRIX<T>& m1, size_t tlr, size_t tlc,
		   size_t brr, size_t brc)
{
	size_t n1=brr-tlr+1, n2=brc-tlc+1;
	m1.resize(n1,n2);
	for (size_t i=0; i<n1; i++)
	{
		size_t p=i+tlr;
		for (size_t j=0; j<n2; j++)
		{
			size_t q=j+tlc;
			m1(i,j)=m(p,q);
		}
	}
}

// split m into m1 and m2
template <class T> 
void split(const MATRIX<T>& m, 
		   MATRIX<T>& m1, MATRIX<T>& m2,
		   size_t m1_count, TFJRowCol rc, bool joint_copy=false)
{
	if (m1_count==0)
	{ m1.resize(0,0); m2=m; return; }

	size_t p=m.size1()-1, q=m.size2()-1,
			pq=m1_count-joint_copy;
	if (rc==rcRow) {
		fjlib::slice_matrix(m,m1,0,0,m1_count-1,q);
		fjlib::slice_matrix(m,m2,pq,0,p,q);
	} else {
		fjlib::slice_matrix(m,m1,0,0,p,m1_count-1);
		fjlib::slice_matrix(m,m2,0,pq,p,q);
	}
}

// copy data m((tlr,tlc)-(brr,brc)) to (dr,dc)
template <class T>
void copy(MATRIX<T>& m, size_t tlr, size_t tlc,
		   size_t brr, size_t brc, size_t dr, size_t dc)
{
	for (size_t i=tlr; i<=brr; i++)
	{
		size_t p=i-tlr+dr;
		for (size_t j=tlc; j<=brc; j++)
		{
			size_t q=j-tlc+dc;
			m(p,q)=m(i,j);
		}
	}
}

// copy data from m1((tlr,tlc)-(brr,brc)) to m(dr,dc)
template <class T>
void copy(MATRIX<T>& m,
		  const MATRIX<T>& m1,
		  size_t tlr, size_t tlc,
		   size_t brr, size_t brc, size_t dr, size_t dc)
{
	for (size_t i=tlr; i<=brr; i++)
	{
		size_t p=i-tlr+dr;
		for (size_t j=tlc; j<=brc; j++)
		{
			size_t q=j-tlc+dc;
			m(p,q)=m1(i,j);
		}
	}
}

// boost resize(m,n,true) doesn't work as expected
// here is corrected version, needs rework
template <class T>
void preserve_resize(MATRIX<T>& m,
				size_t p, size_t q)
{
	size_t op=m.size1(), oq=m.size2();
	MATRIX<T> n(p,q);	
	if (op>p) op=p;
	if (oq>q) oq=q;
	for (size_t i=0; i<op; i++)
		for (size_t j=0; j<oq; j++)
			n(i,j)=m(i,j);
	m=n;
}

// erase m(p0-p1)
template <class T>
void erase(MATRIX<T>& m, 
			size_t p0, size_t p1, TFJRowCol rc)
{
	size_t p=m.size1(), q=m.size2(),
			n=p1-p0+1;
	if (rc==rcRow) {
		if (p1<p-1)
			copy(m,p1+1,0,p-1,q-1,p0,0);
		preserve_resize(m,p-n,q);
	} else {
		if (p1<q-1)
			copy(m,0,p1+1,p-1,q-1,0,p0);
		preserve_resize(m,p,q-n);
	}
}

template <class T>
void split(MATRIX<T>& m, 
		   MATRIX<T>& m2,
		   size_t m_count, TFJRowCol rc, bool joint_copy=false)
{
	if (m_count==0)
	{ m2=m; m.resize(0,0); return; }

	size_t p=m.size1()-1, q=m.size2()-1,
			pq=m_count-joint_copy;
	if (rc==rcRow) {
		fjlib::slice_matrix(m,m2,pq,0,p,q);
		erase(m,m_count,p,rcRow);
	} else {
		fjlib::slice_matrix(m,m2,0,pq,p,q);
		// manually resize m, do it, come on
/*
		MATRIX<T> tmp;
		tmp.resize(p,m_count);
		for (size_t i=0; i<p; i++)
			for (size_t j=0; j<m_count; j++)
				tmp(i,j)=m(i,j);
		m=tmp;
//		m.resize(p,m_count,true);
*/
		erase(m,m_count,q,rcCol);
	}
}

// merge m1, m2 into m
template <class T>
void merge(const MATRIX<T>& m1, 
		   const MATRIX<T>& m2,
		   MATRIX<T>& m,
		   TFJRowCol rc)
{
	if (rc==rcRow) {
		if (m1.size2()!=m2.size2()) throw "can't merge";
		preserve_resize(m,m1.size1()+m2.size1(),m1.size2());
		copy(m,m1,0,0,m1.size1()-1,m1.size2()-1,0,0);
		copy(m,m2,0,0,m2.size1()-1,m2.size2()-1,m1.size1(),0);
	} else {
		if (m1.size1()!=m2.size1()) throw "can't merge";
		preserve_resize(m,m1.size1(),m1.size2()+m2.size2());
		copy(m,m1,0,0,m1.size1()-1,m1.size2()-1,0,0);
		copy(m,m2,0,0,m2.size1()-1,m2.size2()-1,0,m1.size2());
	}
}

// merge m1 into m
template <class T>
void merge(MATRIX<T>& m,
		   const MATRIX<T>& m1,
		   TFJRowCol rc)
{
	if (rc==rcRow) {
		if (m1.size2()!=m.size2()) throw "can't merge";
		int sm=m.size1();
		preserve_resize(m,m1.size1()+sm,m1.size2());
		copy(m,m1,0,0,m1.size1()-1,m1.size2()-1,sm,0);
	} else {
		if (m1.size1()!=m.size1()) throw "can't merge";
		int sm=m.size2();
		preserve_resize(m,m1.size1(),m1.size2()+sm);
		copy(m,m1,0,0,m1.size1()-1,m1.size2()-1,0,sm);
	}
}



} // end of namespace

#endif


/*
	va_list ap;
	va_start(ap,v);
	for (size_t i=0; i<n; i++)
	{
		vec[i]=v;
		v = va_arg(ap, T);
	}
	va_end(ap);	
*/
