//---------------------------------------------------------------------------

#ifndef FJLIBUtilH
#define FJLIBUtilH

// 1.20	some utility math functions, refined

// sqr, sign, swap,
// between, between2, equal, bound,
// relative, balance, balance2,

// lower_pt,upper_pt,lower_seg,upper_seg

namespace fjlib {

// can't believe they don't have sqr
template<class T> inline
const T sqr(const T &a)
{ 
	return a*a; 
}

// return 1 for positve, 0 for zero, -1 for negative
template<class T> inline
const T sgn(const T &a)
{
	return a >= 0 ? (a > 0 ? 1 : 0) : -1;
}

template <class T> inline
void swap(T& a, T& b)
{
	T t=a; a=b; b=t;
}

// returns true if v is between min and max
// include min and max 
template <class T> inline
bool between(const T& v, const T& min, const T& max)
{
	return !((v-min)*(v-max)>0);
}

// returns true if v is between min and max
// exclude max
template <class T> inline
bool between2(const T& v, const T& min, const T& max)
{
	return (!(v<min)&&(v<max));
}

// equal within eps
template <class T> inline
bool equal(const T& v1, const T& v2, const T& eps)
{
	return between(v1,v2-eps,v2+eps);
}

// make v a value between min and max
// min has to less than max
template <class T> inline
T bound(const T& v, const T& min, const T& max)
{
	if (v<min) return min;
	if (v>max) return max;
	return v;
}

// calc relative postion in range
template <class T> inline
T relative(const T& v, const T& _left, const T& _right)
{
	if (_left==_right) return 0.0;
	else return (v-_left)/(_right-_left);
}

// returns average value between left and right according to rel_pos
template <class T> inline
T balance(const T& _left, const T& _right, const T& rel_pos)
{
	return _left*(1-rel_pos)+_right*rel_pos;
}

// returns average value
template <class T> inline
T balance2(const T& _left_low, const T& _right_low,
           const T& _right_top, const T& _left_top,
           const T& hrel_pos, const T& vrel_pos)
{
	return balance(balance(_left_low,_left_top,vrel_pos),
		             balance(_right_low,_right_top,vrel_pos),
			         hrel_pos);
}

// Returns point index of the lower bound of seg segment 
inline size_t	lower_pt(size_t seg)  { return seg; }
// Returns point index of the upper bound of seg segment 
inline size_t	upper_pt(size_t seg)  { return seg+1; }
inline size_t	lower_seg(size_t pt)  { return pt-1; }
inline size_t	upper_seg(size_t pt)  { return pt; }


}  // end of namespace


//---------------------------------------------------------------------------
#endif
