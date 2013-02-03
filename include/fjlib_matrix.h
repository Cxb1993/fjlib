#ifndef FJLIB_MATRIX_H
#define FJLIB_MATRIX_H

#include <vector>

namespace fjlib {

/*!
// A matrix model 
//
// Created, 8/12/2005
!*/
template <class T>
class TFJCustomMatrix {
protected:
	size_t		m_nRow,m_nCol;
	virtual
	void		_resize()=0;
public:
	TFJCustomMatrix(): m_nRow(0), m_nCol(0) {}

	virtual
	void		resize(size_t rows, size_t cols)=0;

	virtual
	const T&	operator()(size_t rows, size_t cols) const=0;
	virtual
	T&			operator()(size_t rows, size_t cols)=0; 

	inline virtual
	size_t		size1() const { return m_nRow; }
	inline virtual
	size_t		size2() const { return m_nCol; }
};

/*!
// A matrix container class with STL compatibility
//
// Created, 8/11/2005
!*/
template <class T>
class TFJMatrix: public TFJCustomMatrix<T>, public std::vector<T> {
public:
	typedef std::vector<T> vec_type;
private:
	inline
	size_t		_total() { return m_nRow*m_nCol; }
	inline
	size_t		_index(int rows, int cols) const { return rows*m_nCol+cols; }
	void		_resize() { vec_type::resize(_total()); }
public:
	TFJMatrix(): TFJCustomMatrix<T>(), vec_type() {}

	void		resize(size_t rows, size_t cols)
	{ m_nRow=rows; m_nCol=cols; _resize(); }

	const T&	operator()(size_t rows, size_t cols) const
	{ return vec_type::operator[](_index(rows,cols)); }
	T&			operator()(size_t rows, size_t cols) 
	{ return vec_type::operator[](_index(rows,cols)); }

	// you need to know what you are doing when calling
	// push_front(), pop_back()
};

}	// end of namespace






#endif