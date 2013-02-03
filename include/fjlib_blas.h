#ifndef FJLIB_BLAS_H
#define FJLIB_BLAS_H

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace fjlib {

namespace ublas=boost::numeric::ublas;

template <class T>
void inv(const ublas::matrix<T>& input, ublas::matrix<T>& inverse) 
{
	using namespace boost::numeric::ublas;
	matrix<T> A(input);
	permutation_matrix<T> P(A.size1());
	lu_factorize(A);
//	lu_factorize(A,P);
	inverse.assign(ublas::identity_matrix<double>(A.size1()));
	lu_substitute<const matrix<T>,matrix<T> >(A,inverse);
}

template <class T>
void solve_axb(const ublas::matrix<T>& Ainput, const ublas::vector<T>& b,
				ublas::vector<T>& x)
{
	using namespace boost::numeric::ublas;
	matrix<T> A(Ainput);
	permutation_matrix<T> P(A.size1());
	lu_factorize(A,P);
	vector<T> tmp=b;
	lu_substitute(A,P,tmp);
	x=tmp;


/*
	ublas::matrix<T> Ainv(A.size1(),A.size2());
	inv(A,Ainv);
	ublas::axpy_prod(Ainv,b,x,true);
*/
}

template <class T>
void solve_trig(const ublas::matrix<T>&A, const ublas::vector<T>& b,
				ublas::vector<T>& x)
{
	x=ublas::solve(A,b,ublas::lower_tag());
}


}	// end of namespace

#endif

