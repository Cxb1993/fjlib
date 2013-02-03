#ifndef FJLIBScaleUserH
#define FJLIBScaleUserH

#include "fjlib_scale.h"

namespace fjlib {

// mapping from node id to nd value
class TFJScale_User: public TFJScale_Base {
protected:
	TFJFunction1p
				func_n2v,func_v2n;
	void*		obj;
public:
	TFJScale_User(): obj(this) {}
	inline
	void		set_user_scale(TFJFunction1p n2v,
								TFJFunction1p v2n)
	{ func_n2v=n2v; func_v2n=v2n; }
	inline
	void		set_obj(void *ob) { obj=ob; }

	inline
	float_t		value_at(float_t pt) 
	{ return func_n2v(obj,&pt); }
	inline
	float_t		index_at(float_t v) 
	{ return func_v2n(obj,&v); }

	inline
	float_t		v_min() { return value_at(nd_min); }
	inline
	float_t		v_max() { return value_at(nd_max); }

};

/* example of construct a unimesh from user define
class TFJScale_UniMesh: public TFJScale_User {
protected:
	float_t		_min,_delta;
	int			_segs;
public:
	inline
	void		set_v_range(float_t low, float_t dx, int segs)
	{ _min=low; _delta=dx; set_nd_range(segs); }
	inline
	float_t		delta() { return _delta; }
	float_t		v_min() { return _min; }
public:
	static
	float_t		n2v(void *obj, float_t* nd_id)
	{ 
		TFJScale_UniMesh* d=(TFJScale_UniMesh*)obj;
		return (*nd_id)*d->delta()+d->v_min();
	}
	static
	float_t		v2n(void *obj, float_t* v)
	{
		TFJScale_UniMesh* d=(TFJScale_UniMesh*)obj;
		return (*v-d->v_min())/d->delta();
	}
	TFJScale_UniMesh() { func_n2v=&n2v; func_v2n=&v2n; }
};
*/

}

#endif