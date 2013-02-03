#ifndef FJLIB_MATRIX_H
#define FJLIB_MATRIX_H

#include <vector>
#include "fjlib_float.h"

// created 8/12/2005, never used
// a small matrix class
// compatible with STL
// modified 3/19/2006, for bds2006

namespace fjlib {

template <class T>
class TFJMatrix_Base {
protected:
	size_t		_rows,_cols;
public:
	TFJMatrix_Base(): _rows(0), _cols(0) {}
	virtual
	void		resize(size_t __rows, size_t __cols)=0;

	virtual
	const T&	operator()(size_t rows, size_t cols) const=0;
	virtual
	T&			operator()(size_t rows, size_t cols)=0; 

	virtual inline
	size_t		rows() const { return _rows; }
	virtual inline
	size_t		cols() const { return _cols; }
};

// derived from STL vector
// these function shouldn't be called in regular basis,
// push_front(), pop_back(), resize(n)
template <class T>
class TFJMatrix: public TFJMatrix_Base<T>, public std::vector<T> {
public:
	typedef std::vector<T> vec_type;
private:
	inline
	size_t		_total() { return _rows*_cols; }
	/// row major index
	inline
	size_t		_index(size_t row, size_t col) const
	{ return row*_cols+col; }
	void		_resize() { vec_type::resize(_total()); }
public:
	TFJMatrix(): TFJMatrix_Base<T>(), vec_type() {}
	TFJMatrix(size_t rows, size_t cols) { resize(rows,cols); }

	void		resize(size_t rows, size_t cols)
	{ _rows=rows; _cols=cols; _resize(); }

	const T&	operator()(size_t rows, size_t cols) const
	{ return vec_type::operator[](_index(rows,cols)); }
	T&			operator()(size_t rows, size_t cols) 
	{ return vec_type::operator[](_index(rows,cols)); }

	// be compatible with other packages
	inline
	size_t 		size1() { return _rows; }
	inline
	size_t 		size2() { return _cols; }
};

/*
typedef TFJMatrix<size_t> matrix_sz;
typedef TFJMatrix<int> matrix_n;
typedef TFJMatrix<float_t> matrix_f;
*/

}	// end of namespace

#endif
