#include "fjlib_interpolation.h"
#include "fjlib_polyroots.h"

namespace fjlib {

bool TFJInterp_Base::splinx(size_t seg, float_t y, vector_f *v,
							float_t eps)
{
//	v->resize(0);

	// construct a polynomial ready for solve
	TFJPolynomial p=get_seg_poly(seg);
	p[0]-=y;

	// solve it
	TFJPolyRoots pr;
	int n=pr.solve(p);

	// if overlapped, very rare but happens
	if (pr.INF)
	{
		push_back(*v,_scale.seg_lower_value(seg));	// use two ends value
		push_back(*v,_scale.seg_upper_value(seg));
		return true;
	}
	// if no root
	if (n<1) return false;

	// handle normal case here
	vector_f rts;
	pr.real_roots(rts);
	bool add=false;
	if (rts.size()>0)
		for (size_t i=0; i<rts.size(); i++)
//			if (between(rts[i],0.0-1e-16,1.0+1e-16)) 
			if (between(rts[i],0.0,1.0)) 
			{ 
				float_t ss=_scale.seg_balance(seg,rts[i]);
				// check uniqueness using eps tolerance
				bool found=false;
				if (eps!=0)
				{
					if (v->size()>0) 
					for (size_t i=0; i<v->size(); i++)
						if (fabs((*v)[i]-ss)<eps) { found=true; break; }
				}
				if (!found)
				{
					push_back(*v,_scale.seg_balance(seg,rts[i]));
					add=true; 
				}
			}
	return add;
}

float_t TFJInterp_Poly::splint(float_t x)
{
	size_t seg=_scale.locate_seg(x);
	float_t frac=_scale.seg_relative(seg,x);
	return _polys[seg].eval(frac);
}

void TFJInterp_Poly::splint(float_t x, vector_f* v)
{
	size_t seg=_scale.locate_seg(x);
	float_t frac=_scale.seg_relative(seg,x);
	int n=_polys[seg].degree();
	_polys[seg].eval(frac,n,v);
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
float_t TFJInterp_Line::splint(float_t x)
{ 
	size_t seg=_scale.locate_seg(x);
	float_t frac=_scale.seg_relative(seg,x);
	return _series.seg_balance(seg,frac);
}

void TFJInterp_Line::splint(float_t x, vector_f* v)
{
	v->resize(2);
	size_t seg=_scale.locate_seg(x);
	float_t frac=_scale.seg_relative(seg,x);
	(*v)[0]=_series.seg_balance(seg,frac);
	(*v)[1]=_series.seg_len(seg)/_scale.seg_len(seg);
}

TFJPolynomial TFJInterp_Line::get_seg_poly(size_t seg)
{
	TFJPolynomial p;
	p.set_degree(1);
	p[0]=_series.seg_lower_value(seg);
	p[1]=_series.seg_upper_value(seg)-p[0];
	return p;
}

/*
bool TFJInterp_Line::splinx(size_t seg, float_t y, vector_f *v)
{
	if (!between(y,_series.seg_lower_value(seg),
					_series.seg_upper_value(seg))) return false;
	float_t frac=_series.seg_relative(seg,y);
	float_t val=_scale.seg_balance(seg,frac);
	// prevent adding same root from nearby segments
	bool ad=true;
	if (v->size()>0) 
		for (size_t i=0; i<v->size(); i++)
			if ((*v)[i]==val) { ad=false; break; }
	if (ad) push_back(*v,val);
	return true;
}
*/


}	// end of namespace

