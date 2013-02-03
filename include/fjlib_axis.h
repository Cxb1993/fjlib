//---------------------------------------------------------------------------
// Axis template
// Axis is an object with uniform spaced segments
// VAxis is an object with non-uniform spacing

#ifndef FJAxisH
#define FJAxisH
//---------------------------------------------------------------------------

#include "fjlib.h"

namespace fjlib {
template <class T>
class TFJAxis_Base {
public:
        virtual size_t	seg_count()=0;
        virtual size_t  pt_count()=0;

        /* Returns point index of the lower bound of seg segment */
        inline size_t   seg_lower_pt(size_t seg) const
        { return seg; }
        /* Returns point index of the upper bound of seg segment */
        inline size_t   seg_upper_pt(size_t seg) const
        { return seg+1; }
        /* Get segment index with the specified point as its lower bound */
        inline size_t   pt_lower_seg(size_t pt) const
        { return pt; }
        /* Get segment index with the specified point as its upper bound */
        inline size_t   pt_upper_seg(size_t pt) const
        { return pt-1; }

		virtual const T pt_pos(size_t pt)=0;
		virtual const T	min0() { return pt_pos(0); }
		virtual const T max0() { return pt_pos(seg_count()); }
        const T  seg_lower(size_t seg)
        { return pt_pos(seg_lower_pt(seg)); }
        inline const T  seg_upper(size_t seg)
        { return pt_pos(seg_upper_pt(seg)); }
        inline const T  seg_center(size_t seg)
        { return (seg_lower(seg)+seg_upper(seg))/2; }
        virtual const T seg_len(size_t seg)
        { return seg_upper(seg)-seg_lower(seg); }

        // native axis marker position manupulation
        inline T        seg_value(size_t seg, const T& rel_v)
        { return seg_len(seg)*rel_v+seg_lower(seg); } // can use balace too
        inline T        seg_rel_value(size_t seg, const T& v)
        { return (v-seg_lower(seg))/seg_len(seg); }

        inline bool_t          inside(const T& v)
        { return between(v,min0(),max0()); }
        bool_t          inside_seg(const T& v, size_t s0, size_t s1)
        {
          if (s1<s0)
            return between(v,seg_lower(s1), seg_upper(s0));
          else
            return between(v,seg_lower(s0), seg_upper(s1));
        }
        inline bool_t	inside_seg(const T& v, size_t seg)
        { return inside_seg(v,seg,seg); }

        virtual int		locate_seg(const T& v)=0;
        virtual size_t  locate_seg_within(const T& v)
        {
          int seg=locate_seg(v);
          if (seg<0) return 0;
          int c=seg_count()-1;
          if (seg>c) return c;
          return seg;
        }
        virtual bool_t	equal_pt(const T& v, size_t pt)
        {
          size_t p=locate_seg(v);
          if (pt_pos(p)==v) return true;
          ++p;
          if ((pt_pos(p)==v) && (p<pt_count())) return true;
          return false;
        }

		virtual void	shift(const T& v)=0;
};

template <class T>
class TFJAxis: public TFJAxis_Base<T> {
protected:
		T				_min, _max;
        size_t          _seg;                   // segment count
        void            _check_seg_count(size_t seg)
        {
          #ifndef _FJLIB_NO_THROW
            if (seg==0) throw std::length_error("can not be zero");
          #endif
        }
public:
        typedef T value_type;
        TFJAxis(): _min(0), _max(1), _seg(1) {}
        ~TFJAxis() {}
        TFJAxis(const T& min_, const T& max_, size_t seg):
                _min(min_), _max(max_), _seg(seg)
        {
          _check_seg_count(seg);
          if (_min>_max) std::swap(_min,_max);
        }
        TFJAxis(const TFJAxis& v): _min(v._min), _max(v._max),
                _seg(v._seg) {}
        TFJAxis&        operator=(const TFJAxis& v)
        {
          if (this!=&v) {
            _min=v._min; _max=v._max; _seg=v._seg; }
          return *this;
        }

        void            set(const T& min_, const T& max_, size_t seg)
        {
          _check_seg_count(seg);
          _min=min_; _max=max_; _seg=seg;
          if (_min>_max) std::swap(_min,_max);
        }
		void			set_seg(size_t seg)
		{
		  _check_seg_count(seg);
		  _seg=seg;
		}

		const T			min0() { return _min; }
		const T			max0() { return _max; }
        size_t			seg_count() { return _seg;}
		size_t			pt_count()  { return _seg+1; }
        inline const T  inc() const { return (_max-_min)/_seg; }

        const T			pt_pos(size_t pt)
        { return _min+pt*inc(); }
        const T			seg_len(size_t seg)
        { return inc(); }

//        T               seg_rel_value(const T&v) const
//        { return seg_rel_value(locate_seg_within(v),v); }

        // linear data interpolation, data type V
        template <class V> inline 
        V               seg_interp_rel(const T& rel,
                                       const V& v0, const V& v1)
        { return v1*rel+v0*(1-rel); }

/* might not be useful in pratical
        template <class V>
        V               seg_interp(const T& v,
                                   const V& v0, const V& v1)
        {
          T rel=seg_rel_value(locate_seg_within(v),v);
          return seg_interp_rel(rel,v0,v1);
        }
*/

        int             locate_seg(const T& v)
        { return std::floor((v-_min)/inc()); }

        // without changing inc()
        bool_t          extend(const T& v, int pin_pos=0);
		void			shift(const T& v) {	_min+=v; _max+=v; }
};

template <class T>
bool_t TFJAxis<T>::extend(const T& v, int pin_pos)
{
  T tmp=_seg+v;
  if (tmp<0) return false;
  switch (pin_pos) {
    case -1:
      _max+=inc()*v;
      break;
    case 1:
      _min-=inc()*v;
      break;
    case 0:
      T mm=_max+inc()*v/2;
      _min-=inc()*v/2;
      _max=mm;
      break;
  }
  _seg=tmp;
  return true;
}

template <class T>
class TFJVAxis: public TFJAxis_Base<T> {
protected:
		std::vector<T>	_vec;
public:
        typedef T				value_type;
		typedef std::vector<T>	vec_type;
        TFJVAxis(): _vec(0) {}
        ~TFJVAxis() {}
		TFJVAxis(const vec_type &vec_) { set(vec_); }
		TFJVAxis(const T& min_, const T& max_, size_t seg) { set_u(min_,max_,seg); }
		TFJVAxis(const TFJVAxis& v) { if (this!=&v) _vec=v._vec; }
        TFJVAxis&        operator=(const TFJVAxis& v)
        {
			if (this!=&v) _vec=v._vec;
			return *this;
        }

		vec_type&		_vector() { return _vec; }
		void			set(const vec_type &vec_) 
		{ if (vec_.size()>0) _vec=vec_; }
        void            set_u(bool is_append, 
							const T& min_, const T& max_, size_t seg)
        {
			size_t vn=_vec.size();
			if (seg<1) return;
			T d=(max_-min_)/seg;			
			size_t mod;
			if (!is_append) { vn=0; mod=0; }
			else mod=1;
			if (vn==0) mod=0;
			_vec.resize(vn+seg+1-mod);		// keep the first node same as the previous one
			if (mod==0) _vec[vn]=min_; 
			for (size_t i=vn+1-mod; i<vn+seg-mod; i++)
				_vec[i]=_vec[i-1]+d;
			_vec[vn+seg-mod]=max_;
        }
		void			set_v1(bool is_append, const T& min_, const T& max_, 
								size_t seg0, size_t seg1)
		{
			size_t vn=_vec.size();
			float_t xmm=max_-min_;
			size_t seg=size_t(seg0*(1.0-1.0/seg1)+1.0+0.5);
			if (seg<3) throw std::length_error("can less than 3");
			size_t mod;
			if (!is_append) { vn=0; mod=0; }
			else mod=1;
			if (vn==0) mod=0;
			_vec.resize(vn+seg+1-mod);
			float_t dx=xmm*(1.0-(float_t)seg/seg0)*2/seg/(seg-1);
			float_t ddx=xmm/seg0;
			if (mod==0) _vec[vn]=min_; 
			_vec[vn+1-mod]=min_+ddx;
			for (size_t i=vn+2-mod; i<vn+seg-mod; i++)
			{
				ddx+=dx;
				_vec[i]=_vec[i-1]+ddx;
			}
			float_t ddd=_vec[vn+seg-1-mod];
			_vec[vn+seg-mod]=max_;
		}

		size_t			pt_count() { return _vec.size(); }
		size_t			seg_count() { return _vec.size()-1; }

		const T			pt_pos(size_t pt) { return _vec[pt]; }

		int				locate_seg(const T& v)
		{
			if (v<min0()) return -1;
			vec_type::iterator p=std::lower_bound(_vec.begin(),_vec.end(),v);
			return pt_upper_seg(p-_vec.begin());
		}

		void			shift(const T& v)
		{  
			using namespace std;
			transform(_vec.begin(),_vec.end(),_vec.begin(),
                bind2nd(plus<T>(),v));
		}
};

typedef TFJAxis<float_t> axis_f;
typedef TFJAxis<int> axis_i;

typedef vector_p<axis_f,2> vector_p2axisf;
typedef vector_p<axis_f,3> vector_p3axisf;

typedef TFJVAxis<float_t> vaxis_f;
typedef TFJVAxis<int> vaxis_i;

typedef vector_p<vaxis_f,2> vector_p2vaxisf;
typedef vector_p<vaxis_f,3> vector_p3vaxisf;

//typedef vector_p2vaxisf grid2_f;


}

#endif

