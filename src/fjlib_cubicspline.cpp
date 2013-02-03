#include "fjlib_cubicspline.h"
#include "fjlib_solver.h"
#include "fjlib_polyroots.h"
 
namespace fjlib {

void TFJCubicSpline::spline()
{
	// prevent same nodes
	if (std::adjacent_find((*_x).begin(),(*_x).end())!=(*_x).end())
		throw "\ncan't have zero length segment from spline() in fjlib_cubicspline.cpp\n";
	// initialize the sparse matrix
	size_t cn=_scale.pt_count();
	ap.resize(cn);	sr.resize(cn);
	al.resize(cn);	ar.resize(cn);
	_y2.resize(cn);
	std::fill(sr.begin(),sr.end(),0.0);
	std::fill(al.begin(),al.end(),0.0);
	std::fill(ar.begin(),ar.end(),0.0);
	// fill out the coefficient except the ends
	for (size_t i=1; i<cn-1; i++)
	{
		float_t dsi=(*_x)[i+1]-(*_x)[i],
				dsi1=(*_x)[i]-(*_x)[i-1];
		al[i]=dsi1/6.0;		
		ar[i]=dsi/6.0;
		ap[i]=(al[i]+ar[i])*2.0;
		sr[i]=((*_y)[i+1]-(*_y)[i])/dsi-
				((*_y)[i]-(*_y)[i-1])/dsi1;
	}
	// special treatment for closed curve, y0=yn-1;
	if ((bt_b==csbtClosed) || (bt_e==csbtClosed))
	{
		size_t dn=cn-1;		// order reduced by 1, careful here!
		float_t ds0=(*_x)[1]-(*_x)[0],
				ds2=(*_x)[cn-1]-(*_x)[cn-2],
				ds3=(*_x)[cn-2]-(*_x)[cn-3];
		ar[0]=ds0/6.0;	al[0]=ds2/6.0;
        ap[0]=(ar[0]+al[0])*2.0;
		sr[0]=((*_y)[1]-(*_y)[0])/ds0-
				((*_y)[cn-1]-(*_y)[cn-2])/ds2;
		_y2.resize(dn);
		TFJCTriDag1D td;
		td.set_data(&ap,&al,&ar,&sr,&_y2);
//		td.al=&al; td.ar=&ar; td.AP=&ap; td.SR=&sr; td.FI=&_y2;
		td.solve();
		_y2.resize(cn);
		_y2[cn-1]=_y2[dn-1];
		return;
	}
	// fill out start b.c.
	if (bt_b==csbtNatural)
	{
		ap[0]=1.0;	sr[0]=yp_b;
	}
	if (bt_b==csbtQuadratic)
	{
		ap[0]=1.0;	ar[0]=-1.0;
	}
	if ((bt_b==csbtClamped) || (bt_b==csbtBessel))
	{
		float_t ds0=(*_x)[1]-(*_x)[0];
		ap[0]=ds0/3.0;	ar[0]=ds0/6.0;
		sr[0]=((*_y)[1]-(*_y)[0])/ds0;
		if (bt_b==csbtBessel)
		{
			float_t ds1=(*_x)[2]-(*_x)[1];
			yp_b=-((*_y)[2]-(*_y)[1])*ds0/ds1/(ds0+ds1)+
				((*_y)[1]-(*_y)[0])*(2*ds0+ds1)/ds0/(ds0+ds1);
		}
		sr[0]-=yp_b;
	}
	// fill out end b.c.
	if (bt_e==csbtNatural)
	{
		ap[cn-1]=1.0;	sr[cn-1]=yp_e;
	}
	if (bt_e==csbtQuadratic)
	{
		ap[cn-1]=1.0;	al[cn-2]=-1.0;
	}
	if ((bt_e==csbtClamped) || (bt_e==csbtBessel))
	{
		float_t ds2=(*_x)[cn-1]-(*_x)[cn-2];
		ap[cn-1]=ds2/3.0;	al[cn-1]=ds2/6.0;
		sr[cn-1]=-((*_y)[cn-1]-(*_y)[cn-2])/ds2;
		if (bt_e==csbtBessel)
		{
			float_t ds3=(*_x)[cn-2]-(*_x)[cn-3];
			yp_e=((*_y)[cn-1]-(*_y)[cn-2])*(ds3+2*ds2)/ds2/(ds2+ds3)-
				((*_y)[cn-2]-(*_y)[cn-3])*ds2/ds3/(ds2+ds3);
		}
		sr[cn-1]+=yp_e;
	}
	TFJTriDag1D td;
	td.set_data(&ap,&al,&ar,&sr,&_y2);
//	td.AL=&al; td.AR=&ar; td.AP=&ap; td.SR=&sr; td.FI=&_y2;
	td.solve();
}

float_t TFJCubicSpline::splint(float_t x)
{
	size_t seg=_scale.locate_seg(x);
	float_t sb,sb1;
	sb=_scale.seg_relative(seg,x);
	sb1=1-sb;
	size_t pt=lower_pt(seg);
	float_t ds=_scale.seg_len(seg);
	return sb1*(*_y)[pt]+sb*(*_y)[pt+1]+
			((sb1*sb1*sb1-sb1)*_y2[pt]+
			(sb*sb*sb-sb)*_y2[pt+1])/6.0*ds*ds;
}

void TFJCubicSpline::splint(float_t x, vector_f* v)
{
	size_t seg=_scale.locate_seg(x);
	float_t sb,sb1;
	sb=_scale.seg_relative(seg,x);
	sb1=1-sb;
	size_t pt=lower_pt(seg);
	float_t ds=_scale.seg_len(seg);
	v->resize(4);
	(*v)[0]=sb1*(*_y)[pt]+sb*(*_y)[pt+1]+
			((sb1*sb1*sb1-sb1)*_y2[pt]+
			(sb*sb*sb-sb)*_y2[pt+1])/6.0*ds*ds;
	(*v)[1]=((*_y)[pt+1]-(*_y)[pt])/ds+
			((1.0-3.0*sb1*sb1)*_y2[pt]+
			(3.0*sb*sb-1.0)*_y2[pt+1])/6.0*ds;
	(*v)[2]=sb1*_y2[pt]+sb*_y2[pt+1];
	(*v)[3]=(_y2[pt+1]-_y2[pt])/ds;
}

TFJPolynomial TFJCubicSpline::get_seg_poly(size_t seg)
{
	float_t ds=_scale.seg_len(seg);;
	float_t ds2=ds*ds;
	TFJPolynomial p;
	p.set_degree(3);

	p[0]=(*_y)[seg];
	p[1]=(*_y)[seg+1]-(*_y)[seg]-(2*_y2[seg]+_y2[seg+1])*ds2/6;
	p[2]=_y2[seg]*ds2/2;
	p[3]=(_y2[seg+1]-_y2[seg])*ds2/6;
	return p;
}

/*
bool TFJCubicSpline::splinx(size_t seg, float_t y, vector_f *v)
{
	TFJPolynomial p=get_seg_poly(seg);
	p[0]-=y;
	TFJPolyRoots pr;
	pr.solve(p);
	vector_f rts;
	pr.real_roots(rts);
	bool add=false;
	if (rts.size()>0)
		for (size_t i=0; i<rts.size(); i++)
			if (between(rts[i],0.0,1.0)) 
			{ 
				push_back(*v,_scale.seg_balance(seg,rts[i]));
				add=true; 
			}
	return add;
}
*/

}	// end of namespace

