#ifndef TFJAppStaggeredAxisH
#define TFJAppStaggeredAxisH

#include "fjlib_scale.h"
#include <vector>

namespace fjlib {

// 5.25 Derived new scale class to have the uniform
//		structure to be feed into staggered template class

// staggered template interface
class IFJScale_Staggered_Template {
public:
	// make the scale with double refined mesh
	virtual void	refine()=0;
	// spring the scale accoording to the style
	// style=0, half indent
	// style=1, full indent
	virtual void	spring(int style)=0;
};

class TFJScale_Staggered_UniMesh: public TFJScale_UniMesh,
								public IFJScale_Staggered_Template {
public:
	void	refine()
	{
		_delta/=2.0; _segs*=2;
	}
	void	spring(int style)
	{
		if (style==0)
		{
			_min-=delta()/2.0; _segs+=1;
		}
		else if (style==1)
		{
			_min-=delta(); _segs+=2;
		}
	}
};

class TFJScale_Staggered_Mesh: public TFJScale_Mesh,
								public IFJScale_Staggered_Template {
private:
	void			make_vec(vector_f& target, int style)
	{
		int nseg=seg_count();
		float_t lenL=seg_len(0),
				lenR=seg_len(nseg-1);
		// make half indent left
		if (style==0)
		{
			target.resize(nseg+2);
			target[0]=min()-lenL/2.0;
			for (int i=0; i<nseg; i++)
				target[i+1]=seg_center_value(i);
			target[nseg+1]=max()+lenR/2.0;
		}
		// make full indent left
		else if (style==1)
		{
			target.resize(nseg+3);
			target[0]=min()-lenL;
			for (int i=0; i<=nseg; i++)
				target[i+1]=seg_lower_value(i);
			target[nseg+2]=max()+lenR;
		}
	}
public:
	void	refine()
	{
		int sold=seg_count();
		if (sold<1) return;
		int snew=sold*2;
		vector_f b(snew+1);
		b[0]=(*vec)[0];
		for (size_t i=0; i<sold; i++)
		{
			b[i*2+1]=seg_center_value(i);
			b[i*2+2]=seg_upper_value(i);
		}
		*vec=b;
	}
	void	spring(int style)
	{
		vector_f v;
		make_vec(v,style);
		values=v;
	}
};

typedef TFJGroup_XYZ<TFJScale_Staggered_UniMesh> TFJAxis_Staggered_UniMesh;
typedef TFJGroup_XYZ<TFJScale_Staggered_Mesh> TFJAxis_Staggered_Mesh;

/*!
//	2D Staggered structured mesh position type
*/
enum TFJStaggeredPosType {
	sptO=0,sptU,sptV,sptP,sptN
};

/*!
//	2D Staggered structured mesh helper class
//
//	5.25 NonUniform mesh support added 
//	The pain is that the treatment between uni and non-uniform
//	mesh is so different that they have to splited into 
//	different class instead of still using template
*/
template <class ST=TFJScale_Staggered_UniMesh>
class TFJStaggeredAxis_Helper {
public:
	typedef ST						mesh_type;
	typedef TFJGroup_XYZ<mesh_type>	axis_type;
	typedef TFJStaggeredPosType		staggered_type;
	typedef std::vector<axis_type>	axis_vec_type;
protected:
	axis_vec_type	vaxis;
	/// Init internal setup
	virtual
	void			setup(const axis_type& _axis)
	{
		vaxis.resize(5);
		for (size_t i=0; i<5; i++) vaxis[i].redim(2);

		mesh_type xh,xf,yh,yf;
		xh=_axis.x(); xh.spring(0);
		xf=_axis.x(); xf.spring(1);
		yh=_axis.y(); yh.spring(0);
		yf=_axis.y(); yf.spring(1);

		// sptO, copy
		vaxis[sptO].x()=_axis.x();
		vaxis[sptO].y()=_axis.y();
		// sptU
		vaxis[sptU].x()=xh;
		vaxis[sptU].y()=yf;
		// sptV
		vaxis[sptV].x()=xf;
		vaxis[sptV].y()=yh;
		// sptP
		vaxis[sptP].x()=xf;
		vaxis[sptP].y()=yf;
		// sptN
		vaxis[sptN].x()=xh;
		vaxis[sptN].y()=yh;
	}
public:
	TFJStaggeredAxis_Helper(): _curStag(sptO) {}
	TFJStaggeredAxis_Helper(const axis_type& ax): _curStag(sptO)
	{ setup(ax); }
	void		set(const axis_type& ax) { setup(ax); }

	inline
	axis_type&		axis(staggered_type stag)
	{ return vaxis[stag]; }

	inline const
	axis_type&		axis(staggered_type stag) const
	{ return vaxis[stag]; }

	inline const
	size_t			xseg_count(staggered_type stag) const
	{ return vaxis[stag].x().seg_count(); }

	inline const
	size_t			yseg_count(staggered_type stag) const
	{ return vaxis[stag].y().seg_count(); }

	inline const 
	float_t			xseg_len(staggered_type stag, int m) const 
	{ return vaxis[stag].x().seg_len(m); }

	inline const 
	float_t			yseg_len(staggered_type stag, int n) const 
	{ return vaxis[stag].y().seg_len(n); }

	inline const
	float_t			dx(staggered_type stag, int m) const
	{ return xseg_len(stag,m); }
	inline const
	float_t			dy(staggered_type stag, int n) const
	{ return yseg_len(stag,n); }

	inline const
	float_t			xmin(staggered_type stag) const
	{ return vaxis[stag].x().min(); }
	inline const
	float_t			ymin(staggered_type stag) const
	{ return vaxis[stag].y().min(); }
	inline const
	float_t			xmax(staggered_type stag) const
	{ return vaxis[stag].x().max(); }
	inline const
	float_t			ymax(staggered_type stag) const
	{ return vaxis[stag].y().max(); }

	inline const
	float_t			xseg_center(staggered_type stag, int m) const
	{ return vaxis[stag].x().seg_center_value(m); }
	inline const
	float_t			yseg_center(staggered_type stag, int n) const
	{ return vaxis[stag].y().seg_center_value(n); }
	inline const
	float_t			xp(staggered_type stag, int m) const
	{ return xseg_center(stag,m); }
	inline const
	float_t			yp(staggered_type stag, int n) const
	{ return yseg_center(stag,n); }
protected:
	staggered_type	_curStag;
public:	// make it easier by setting current stag type
	void			set_Stag(staggered_type stag)
	{ _curStag=stag; }
	inline const
	size_t			xc() const
	{ return xseg_count(_curStag); }
	inline const
	size_t			yc() const
	{ return yseg_count(_curStag); }
	inline const
	float_t			dx(int m) const
	{ return xseg_len(_curStag,m); }
	inline const
	float_t			dy(int n) const
	{ return yseg_len(_curStag,n); }
	inline const
	float_t			xmin() const
	{ return xmin(_curStag); }
	inline const
	float_t			xmax() const
	{ return xmax(_curStag); }
	inline const
	float_t			ymin() const
	{ return ymin(_curStag); }
	inline const
	float_t			ymax() const
	{ return ymax(_curStag); }
	inline const
	float_t			xp(int m) const
	{ return xseg_center(_curStag,m); }
	inline const
	float_t			yp(int n) const
	{ return yseg_center(_curStag,n); }
};

typedef TFJStaggeredAxis_Helper<TFJScale_Staggered_UniMesh> 
						TFJStagAxisHelper_UniMesh;
typedef TFJStaggeredAxis_Helper<TFJScale_Staggered_Mesh> 
						TFJStagAxisHelper_Mesh;

/*
class TFJStaggered_UniMesh2D: 
		public TFJStaggered_Base2D<TFJScale_UniMesh> {
protected:
	void			setup(const axis_type& _axis) 
	{
		TFJStaggered_Base2D<TFJScale_UniMesh>::setup(_axis);
		float_t xmin=_axis.x().min(),
				ymin=_axis.y().min();
		size_t	xseg=_axis.x().seg_count(),
				yseg=_axis.y().seg_count();
		float_t dx=_axis.x().delta(),
				dy=_axis.y().delta();
		float_t dx2=dx/2.0,
				dy2=dy/2.0;

		vaxis[sptO].x().set_range(xmin,dx,xseg);
		vaxis[sptO].y().set_range(ymin,dy,yseg);
		vaxis[sptU].x().set_range(xmin-dx2,dx,xseg+1);
		vaxis[sptU].y().set_range(ymin-dy,dy,yseg+2);
		vaxis[sptV].x().set_range(xmin-dx,dx,xseg+2);
		vaxis[sptV].y().set_range(ymin-dy2,dy,yseg+1);
		vaxis[sptP].x().set_range(xmin-dx,dx,xseg+2);
		vaxis[sptP].y().set_range(ymin-dy,dy,yseg+2);
		vaxis[sptN].x().set_range(xmin-dx2,dx,xseg+1);
		vaxis[sptN].y().set_range(ymin-dy2,dy,yseg+1);
	}
public:
	TFJStaggered_UniMesh2D(const axis_type& ax) { setup(ax); }
	inline const
	float_t			dx() const
	{ return vaxis[0].x().delta(); }
	inline const
	float_t			dy() const
	{ return vaxis[0].y().delta(); }
};

class TFJStaggered_Mesh2D: 
		public TFJStaggered_Base2D<TFJScale_Mesh> {
private:
	void			make_vec(const mesh_type& src,
							vector_f& target, int style)
	{
		int nseg=src.seg_count();
		float_t lenL=src.seg_len(0),
				lenR=src.seg_len(nseg-1);
		// make half indent left
		if (style==0)
		{
			target.resize(nseg+2);
			target[0]=src.min()-lenL/2.0;
			for (int i=0; i<nseg; i++)
				target[i+1]=src.seg_center_value(i);
			target[nseg+1]=src.max()+lenR/2.0;
		}
		// make full indent left
		else if (style==1)
		{
			target.resize(nseg+3);
			target[0]=src.min()-lenL;
			for (int i=0; i<=nseg; i++)
				target[i+1]=src.seg_lower_value(i);
			target[nseg+2]=src.max()+lenR;
		}
	}
protected:
	void			setup(const axis_type& _axis) 
	{
		TFJStaggered_Base2D<TFJScale_Mesh>::setup(_axis);

		float_t xmin=_axis.x().min(),
				ymin=_axis.y().min();
		size_t	xseg=_axis.x().seg_count(),
				yseg=_axis.y().seg_count();

		vector_f ox=*_axis.x().get_vector(),
				 oy=*_axis.y().get_vector();
		vector_f txf,tyf,txh,tyh;
		make_vec(_axis.x(),txh,0);
		make_vec(_axis.x(),txf,1);
		make_vec(_axis.y(),tyh,0);
		make_vec(_axis.y(),tyf,1);

		// sptO, copy
		vaxis[sptO].x().set_values(ox);
		vaxis[sptO].y().set_values(oy);
		// sptU
		vaxis[sptU].x().set_values(txh);
		vaxis[sptU].y().set_values(tyf);
		// sptV
		vaxis[sptV].x().set_values(txf);
		vaxis[sptV].y().set_values(tyh);
		// sptP
		vaxis[sptP].x().set_values(txf);
		vaxis[sptP].y().set_values(tyf);
		// sptN
		vaxis[sptN].x().set_values(txh);
		vaxis[sptN].y().set_values(tyh);
	}
public:
	TFJStaggered_Mesh2D(const axis_type& ax) { setup(ax); }
};
*/
/*
class TFJStaggered_Mesh2D {
public:
	typedef TFJScale_Mesh			mesh_type;
	typedef TFJGroup_XYZ<mesh_type>	axis_type;
	typedef TFJStaggeredPosType		staggered_type;
	typedef std::vector<axis_type>	axis_vec_type;
protected:
	axis_vec_type	vaxis;
	void			setup(const axis_type& _axis) 
	{
		float_t xmin=_axis.x().min(),
				ymin=_axis.y().min();
		size_t	xseg=_axis.x().seg_count(),
				yseg=_axis.y().seg_count();
		float_t dx=_axis.x().delta(),
				dy=_axis.y().delta();
		float_t dx2=dx/2.0,
				dy2=dy/2.0;

		vaxis.resize(5);
		for (size_t i=0; i<5; i++) vaxis[i].redim(2);
		vaxis[sptO].x().set_range(xmin,dx,xseg);
		vaxis[sptO].y().set_range(ymin,dy,yseg);
		vaxis[sptU].x().set_range(xmin-dx2,dx,xseg+1);
		vaxis[sptU].y().set_range(ymin-dy,dy,yseg+2);
		vaxis[sptV].x().set_range(xmin-dx,dx,xseg+2);
		vaxis[sptV].y().set_range(ymin-dy2,dy,yseg+1);
		vaxis[sptP].x().set_range(xmin-dx,dx,xseg+2);
		vaxis[sptP].y().set_range(ymin-dy,dy,yseg+2);
		vaxis[sptN].x().set_range(xmin-dx2,dx,xseg+1);
		vaxis[sptN].y().set_range(ymin-dy2,dy,yseg+1);
	}
public:
	TFJStaggered_Mesh2D() {}
	TFJStaggered_Mesh2D(const axis_type& ax) { setup(ax); }
	void		set(const axis_type& ax) { setup(ax); }

	inline
	axis_type&		axis(staggered_type stag)
	{ return vaxis[stag]; }

	inline const
	axis_type&		axis(staggered_type stag) const
	{ return vaxis[stag]; }

	inline const
	size_t			xseg_count(staggered_type stag) const
	{ return vaxis[stag].x().seg_count(); }

	inline const
	size_t			yseg_count(staggered_type stag) const
	{ return vaxis[stag].y().seg_count(); }
	inline const
	float_t			dx() const
	{ return vaxis[0].x().delta(); }
	inline const
	float_t			dy() const
	{ return vaxis[0].y().delta(); }


};
*/
}		// end of namespace

#endif
