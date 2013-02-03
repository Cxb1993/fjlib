#ifndef FJLIB_RANGE_H
#define FJLIB_RANGE_H

#include "fjlib_xyz.h"

namespace fjlib {

// created 3/7/2006

template <class T>
struct range {
	T min,max;
	range() {}
	range(const T& _min, const T& _max) {
		min=_min; max=_max;
	}
	range(const range& _range) {
		if (this!=*_range) {
			min=_range.min; 
			max=_range.max;
		}
	}
	virtual inline
	const T length() const { return max-min; } 
	virtual inline
	const T	center() const { return length()/2+min; }
	/// check if value v is within this range
	virtual inline
	bool	within(const T& v) const 
	{ return ((v>=min) && (v<=max)); }
	/// bound the value v within this range
	virtual
	T		bound(const T& v) const 
	{
		if (v<min) return min;
		if (v>max) return max;
		return v;
	}
	/// expand the range by starting from center to radius of v/2
	virtual 
	void	expand(const T& _length) {
		T dv=_length/2,m=center();
		min=m-dv; max=m+dv;
	}
};

template <class T>
std::ostream& operator<<(std::ostream& out, const range<T>& v)
{
	out << v.min << ' ' << v.max;
	return out;
}

template <class T>
std::istream& operator>>(std::istream& in, const range<T>& v)
{
	in >> v.min >> v.max;
	return in;
}

typedef range<size_t> range_sz;
typedef range<int> range_n;
typedef range<float_t> range_f;

template <class T>
struct range2: public xy<range<T> >  {
	typedef T					value_type;
	typedef range<T>			range_type;
	typedef xy<T>				point_type;
	typedef xy<range<T> > 		range_xy_type;

	range2(): range_xy_type() {}
	range2(const range2& v): range_xy_type(v) {}
	range2(const range_type& _x,
			const range_type& _y): 
			range_xy_type(_x,_y) {}
	range2(const point_type& left_bottom, 
			const point_type& right_top) {
		x=range_type(left_bottom.x,right_top.x);
		y=range_type(left_bottom.y,right_top.y);
	}
	range2(const T& xmin, const T& xmax,
			const T& ymin, const T& ymax) {
		x=range_type(xmin,xmax);
		y=range_type(ymin,ymax);
	}
	const point_type	left_bottom() const {
		return point_type(x.min,y.min);
	}
	const point_type 	right_top() const {
		return point_type(x.max,y.max);
	}
	const point_type 	center() const {
		return point_type(x.center(),y.center());
	}
	virtual inline
	bool				within(const T& _x, const T& _y)
	{ return (x.within(_x) && y.within(_y)); }
	virtual inline
	bool				within(const point_type& p)
	{ return within(p.x,p.y); }
	virtual 		
	point_type			bound(const point_type& p)
	{ point_type v; v.x=x.bound(p.x); v.y=y.bound(p.y); return v; }
};

typedef range2<size_t> range2_sz;
typedef range2<int> range2_n;
typedef range2<float_t> range2_f;

}	// end of namespace

#endif
