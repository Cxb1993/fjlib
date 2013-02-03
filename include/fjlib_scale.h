#ifndef FJLIBScaleH
#define FJLIBScaleH

#include "fjlib.h"
//#include "fjlib_function.h"

namespace fjlib {

/*! 
// Series base class: \n
// A small wrapper on a series of float data
// supplies several functions to easy access
// the segment and nodes value. \n
*/
class TFJSeries_Base {
public:
/*
	TFJSeries_Base() {}
	TFJSeries_Base(const TFJSeries_Base& v) {}
*/
	/// Returns the total number of segments
	virtual 
	size_t		seg_count() const=0;
	/// Returns the total number of nodes, usually is seg_count()+1
	virtual 
	size_t		pt_count() const { return seg_count()+1; }
	/// Returns the scale value at node pt
	virtual 
	float_t		value_at(size_t pt) const=0;
/* helper functions */
	/// Returns the left scale value of segment |-| 
	inline 
	float_t		seg_lower_value(size_t seg) const
	{ return value_at(lower_pt(seg)); }
	/// Returns the right scale value of segment |-| 
	inline 
	float_t		seg_upper_value(size_t seg) const
	{ return value_at(upper_pt(seg)); }
	inline 
	float_t		seg_relative(size_t seg, float_t v) const
	{ return relative(v,seg_lower_value(seg),seg_upper_value(seg)); }
	inline 
	float_t		seg_balance(size_t seg, float_t v) const
	{ return balance(seg_lower_value(seg),seg_upper_value(seg),v); }
	inline 
	float_t		seg_center_value(size_t seg) const
	{ return seg_balance(seg,0.5); }
	/// Returns the length of the segment
	inline virtual
	float_t		seg_len(size_t seg) const
	{ return seg_upper_value(seg)-seg_lower_value(seg); }
};

/*! 
// Series class for vector: 
//
// Wrap a series helper class on a vector \n
// Essentially this class can be treated as a general vector helper
*/ 
class TFJSeries: public TFJSeries_Base {
protected:
	vector_f*	vec;
public:
	TFJSeries(): vec(NULL) {}
	TFJSeries(const TFJSeries& v):TFJSeries_Base(v) {
		if (this!=&v) vec=v.vec; }
	inline
	void		link_vector(vector_f* _vec) { vec=_vec; }
	inline
	size_t		seg_count() const { return vec->size()-1; }
	inline
	float_t		value_at(size_t pt) const { return (*vec)[pt]; }
};

/*!
// Scale base class:
// Helper class which help you to access the data on 
// a 1D scale which is divided by nodes. 
// The index of nodes are fixed to be integals and start 
// from zero but the scale of each node can be accessed
// by calling value_at(). This is defined in its derivative
// classes.
// The Scale is normally monotonically increasing and 
// bounded by its min and max
//
// 1.20	revised, template deleted
// 5.25 refine() added to facilicate staggered mesh
// 1.5	correct warning of base class need to be explicitly
//		called in copy constructor
*/

class TFJScale_Base: public TFJSeries_Base {
/* basic functions */
public:
	TFJScale_Base() {}
	TFJScale_Base(const TFJScale_Base& v):TFJSeries_Base(v) {}
	/// Returns the minimium of the scale
	virtual 
	float_t		min() const { return value_at(0); }
	/// Returns the maximium of the scale
	virtual
	float_t		max() const { return value_at(seg_count()); }
	virtual
	float_t		len() const { return max()-min(); }
/* more functions */
	/// Check if v is inside the whole scale
	inline
	bool		inside(float_t v) const
	{ return between(v,min(),max()); }
	/// Check if v is inside segments [seg1,seg2]
	inline
	bool		inside_segs(float_t v, size_t seg1, size_t seg2) const
	{ return between(v,seg_lower_value(seg1),
						seg_upper_value(seg2)); }
	/// Check if v is inside the segment seg
	inline
	bool		inside_seg(float_t v, size_t seg) const
	{ return inside_segs(v,seg,seg); }
/* advanced functions */
	/// Locate the nearest segment from value v
    virtual		
	size_t		locate_seg(float_t v) const=0;
};

/*!
// Scale class with uniform grid spacing \n
// Call set_range() to define the grids
*/
class TFJScale_UniMesh: public TFJScale_Base {
protected:
	float_t		_min,_delta;
	size_t		_segs;
public:
	/// Default Constructor
	TFJScale_UniMesh(): _min(0),_delta(0),_segs(1) {}
	/// Constructor using set_range()
	TFJScale_UniMesh(float_t low, float_t dx, size_t segs)
	{ set_range(low,dx,segs); }
	/// Define the grids by min, spacing and seg number
	inline
	void		set_range(float_t low, float_t dx, size_t segs)
	{ _min=low; _delta=dx; _segs=segs; }
	inline
	/// Returns the spacing
	float_t		delta() const { return _delta; }
	float_t		min() const { return _min; }
public:
	inline
	float_t		seg_len(size_t seg) const
	{ return delta(); }
	inline
	size_t		seg_count() const { return _segs; }
	inline
	float_t		value_at(size_t pt) const
	{ return _min+pt*_delta; }
	size_t		locate_seg(float_t v) const
	{ 
		if (fabs(_delta)<FJLIB_FLOAT_EPS) 
			throw "\ndelta can't be zero from locate_seg(float_t) in fjlib_scale.h\n";
		int rel=(int)((v-_min)/_delta); 
		return bound(rel,0,(int)seg_count()-1);
	}
public:		
	/// Copy constructor
	TFJScale_UniMesh(const TFJScale_UniMesh& v): TFJScale_Base(v)
	{ 
		if (&v!=this) {
			_min=v._min; _delta=v._delta; _segs=v._segs; }
	}
	/// Operator =, assign
	TFJScale_UniMesh&
				operator=(const TFJScale_UniMesh& v)
	{
		if (this!=&v) {
			_min=v._min; _delta=v._delta; _segs=v._segs; }
        return *this;
	}
};

/*!
// Scale class with non-uniform grid spacing
// with data linked with external vector
// 5.23 GhostMesh and Mesh added
*/
class TFJScale_GhostMesh: public TFJScale_Base {
protected:
	vector_f*	vec;
private:
	void		_add_range(float_t dx, size_t segs)
	{
		if (vec->size()<1) 
			throw "\nneed at least one node in fjlib_scale.h\n";
		size_t on=vec->size();
		size_t nn=on+segs;
		vec->resize(nn);
		float_t low=(*vec)[on-1];
		for (size_t i=0; i<segs; i++)
//			values[on+i]=values[on+i-1]+dx; // bigger error??
			(*vec)[on+i]=low+dx*(i+1);
	}
public:
	TFJScale_GhostMesh(): vec(NULL) {}
	TFJScale_GhostMesh(const TFJScale_GhostMesh& v):TFJScale_Base(v)
	{ 
		if (&v!=this) { vec=v.vec; }
	}
	vector_f*	get_vector() const { return vec; }
	/// Link with an external vector
	virtual
	void		link_vector(vector_f* _vec) { vec=_vec; }
	inline
	size_t		seg_count() const { return vec->size()-1; }

	inline
	float_t		value_at(size_t pt) const { return (*vec)[pt]; }
	/// Values can be copied from a vector
	void		set_values(vector_f _val) { *vec=_val; }
	/// Values can be added by inserting range pieces
	void		set_range(float_t low, float_t dx, size_t segs)
	{
		vec->resize(1);
		(*vec)[0]=low;
		_add_range(dx,segs);
	}
	/// Values can be added by inserting range pieces
	void		add_range(float_t dx, size_t segs)
	{
//		int n=vec->size();
//		if (n<1) throw "\n add_range can't be called before create_range in fjlib_scale.h\n";
		_add_range(dx,segs);
	}
	// Locate valve v
	size_t		locate_seg(float_t v) const
	{ 
		size_t i=1;
		if (v<=min()) return 0;
		while (i<pt_count())
		{
			if (v<=value_at(i)) return lower_seg(i);
			i++;
		}
		return seg_count()-1;
	}
};

/*!
// Scale class with non-uniform grid spacing \n
// value data is store internally
//
*/
class TFJScale_Mesh: public TFJScale_GhostMesh {
protected:
	/// Grid scale value vector
	vector_f	values;
private:
	void		_link_vec()
	{ vec=&values; }
public:
	TFJScale_Mesh() { _link_vec(); }
	/// Override default copy which will point vec to the 
	/// same source vec
	TFJScale_Mesh(const TFJScale_Mesh& v): TFJScale_GhostMesh(v)
	{ 
		if (&v!=this) {
			values=v.values; _link_vec();
		}
	}
	/// Override default assign operator
	/// If not present, TFJGroup_XYZ will fail
	TFJScale_Mesh&
				operator=(const TFJScale_Mesh& v)
	{
		if (&v!=this) {
			values=v.values; _link_vec();
		}
		return *this;
	}

	void		link_vector(vector_f* _vec) 
	{ throw "\nlink_vector is only usable for ghost_mesh in fjlib_scale.h\n"; }
};

}	// end of namespace

#include "fjlib_group.h"

namespace fjlib {

/// Multi-Dimension UniMesh Scale
typedef TFJGroup_XYZ<TFJScale_UniMesh>	TFJAxis_UniMesh;
/// Multi-Dimension Non-UniMesh Scale
typedef TFJGroup_XYZ<TFJScale_Mesh>		TFJAxis_Mesh;
typedef TFJGroup_XYZ<TFJScale_GhostMesh>
										TFJAxis_GMesh;

}	// end of namespace

#endif


