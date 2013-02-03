#include "fjapp_vof_staggered.h"
 
namespace fjlib {

template <class ST,	class CT>
void TFJVOFGen_Staggered<ST,CT>::set_data(
				axis_type* _axis, curve_type* _curve)
{
	axis=_axis; curve=_curve; 
	naxis.redim(2);
	for (size_t i=0; i<2; i++)
	{
		naxis.at(i)=axis->at(i);
		naxis.at(i).refine();
	}
	vof.set_data(&naxis,curve);
}

template <class ST,	class CT>
void TFJVOFGen_Staggered<ST,CT>::vof_gen(bool fill)
{
	vof.vof_gen();
	vof.calc_vofs(mat);
	if (fill)
		vof.fill_vofs(mat);
}

template <class ST,	class CT>
float_t TFJVOFGen_Staggered<ST,CT>::cell_area(size_t m, size_t n)
{
	if (!vof.cylinder)
		return naxis.x().seg_len(m)*naxis.y().seg_len(n);
	else
		return sqr(naxis.x().seg_upper_value(m))-
				sqr(naxis.x().seg_lower_value(m));
/*
	if (!vof.cylinder)
		return naxis.x().delta()*naxis.y().delta();
	else
		return sqr(naxis.x().seg_upper_value(m))-
				sqr(naxis.x().seg_lower_value(m));
*/
}

template <class ST,	class CT>
float_t TFJVOFGen_Staggered<ST,CT>::get_staggered_vof(size_t m, size_t n, 
							staggered_type stag_loc)
{
	float_t sum=0.0,area=0.0;
	int nxseg,nyseg;
	for (size_t cell=0; cell<7; cell+=2)
	{
		staggered_mn(m,n,stag_loc,cell,nxseg,nyseg);
		nxseg=bound(nxseg,0,(int)naxis.x().seg_count()-1);
		nyseg=bound(nyseg,0,(int)naxis.y().seg_count()-1);
		sum+=mat(nxseg,nyseg)*cell_area(nxseg,nyseg);
		area+=cell_area(nxseg,nyseg);
	}
	return sum/area;
}

template <class ST,	class CT>
bool TFJVOFGen_Staggered<ST,CT>::is_colored(size_t m, size_t n,
					staggered_type stag_loc,
					size_t cell_loc)
{
	int nxseg,nyseg;
	staggered_mn(m,n,stag_loc,cell_loc,nxseg,nyseg);
	bool in=between(nxseg,0,(int)naxis.x().seg_count()-1)
			&& between(nyseg,0,(int)naxis.y().seg_count()-1);
	if (!in) return false;
	else return vof.is_colored(nxseg,nyseg);
}

/*
template <class ST, class CT>
void TFJVOFGen_Staggered<ST,CT>::fill_vofs(matrix_f& mat)
{
	size_t m=axis->x().seg_count();
	size_t n=axis->y().seg_count();
	mat.resize(m,n);
	for (size_t i=0; i<m ;i++)
		for (size_t j=0; j<n; j++)
			mat(i,j)=vof.get_staggered_vof(i,j,sptO);
}
*/

template <class ST, class CT>
void TFJVOFGen_Staggered<ST,CT>::convert_xtoseg(const vector_f& v1,
											vector_n & v2)
{
	size_t xn=v1.size();
	v2.resize(xn);
	if (xn<1) return;
	for (size_t i=0; i<xn; i++)
		v2[i]=this->axis->x().locate_seg(v1[i]);	
}

template <class ST,	class CT>
void TFJVOFGen_Staggered<ST,CT>::fill_unknowns(matrix_n& umat)
{
	size_t xn=axis->x().seg_count();
	size_t yn=axis->y().seg_count();
	umat.resize(xn,yn);
	reset_matrix(umat,1);

	size_t j=0;
	vector_f xc;
	vector_n xs; vector_b xb;
	while (j<yn)
	{
		// xc - coord, xs - segment, xb - over middle x
		// make a complete region segments
		vector_f yrx;
		vof.get_y_roots_x(this->axis->y().seg_center_value(j),yrx);
		size_t xnn=yrx.size();
		if (xnn<1) { j++; continue; }
		xc.resize(xnn+2);
		for (size_t i=0; i<xnn; i++)
			xc[i+1]=yrx[i]; 
		xc[0]=this->axis->x().min();
		xc[xnn+1]=this->axis->x().max();
		xnn+=2;
		convert_xtoseg(xc,xs);
		xb.resize(xnn);
		for (size_t i=0; i<xnn; i++)
			xb[i]=(xc[i]>this->axis->x().seg_center_value(xs[i]));
		
		int m,n;
		// find out vacant block known site
		for (size_t i=0; i<xnn; i+=2)
		{
			// don't want to hit the last end
			if (i+1>=xnn) continue;
			
			m=xs[i]; n=xs[i+1];
			if (xb[i]) m++;
			if (!(xb[i+1])) n--; 
			if (n>=m) {
				for (size_t k=m; k<=(size_t)n; k++)
					umat(k,j)=0;		// block knowns
			}	
		}
		j++;
	}
}
	

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
template <class ST,	class CT>
void TFJVOFGen_Staggered_ST<ST,CT>::get_cell_st(size_t m,size_t n, 
								float_t& ust, float_t& vst)
{
	ust=0.0;	vst=0.0;
//	if (!vof.is_colored(m,n)) return 0.0;

	// can use the above code but it's been tagged already
	float_t test=this->mat(m,n);
	if ((this->mat(m,n)==0) || (this->mat(m,n)==1)) return; 
	this->vof.is_colored(m,n);

	TFJPolynomial p,q,p1,p2,q1,q2,pl;
	float_t ds,eps=TFJVOFGen_Staggered<ST,CT>::get_vof().eps;
	vector_f arc;
	typename TFJVOFGen_Staggered<ST,CT>::vof_curve_type&
					vc=this->vof.vof_curve(&arc);
	for (size_t i=0; i<vc.seg_count(); i++)
	{
		ds=vc.arcs().seg_len(i);
		// make sure don't get zero width segment
		if (ds<eps) continue;
		p=vc.interps().x().get_seg_poly(i);

//		if (p.degree()<2) continue;
		if ((arc[i]<0) || (arc[i+1]<0)) continue;

 		q=vc.interps().y().get_seg_poly(i);
		p1=p>>1;	p2=p1>>1;
		q1=q>>1;	q2=q1>>1;

        pl=(p2*p-q1*q1)<<1;
		ust+=pl.eval_diff(0.0,1.0)/ds;
		pl=(p*q*q1)<<1;
		ust-=Bo*pl.eval_diff(0.0,1.0);

        pl=(p*q2+p1*q1)<<1;
		vst+=pl.eval_diff(0.0,1.0)/ds;
		pl=(q*p*p1)<<1;
		vst+=Bo*pl.eval_diff(0.0,1.0);
	}
}

template <class ST,	class CT>
void TFJVOFGen_Staggered_ST<ST,CT>::get_cell_stg(size_t m,size_t n, 
								float_t& ust, float_t& vst)
{
	ust=0.0;	vst=0.0;
//	if (!vof.is_colored(m,n)) return 0.0;

	// can use the above code but it's been tagged already
	float_t test=this->mat(m,n);
	if ((this->mat(m,n)==0) || (this->mat(m,n)==1)) return; 
	this->vof.is_colored(m,n);

	TFJPolynomial p,q,p1,p2,q1,q2,pl,sg,sg1;
	float_t sga,sgb;
	float_t ds,eps=TFJVOFGen_Staggered<ST,CT>::get_vof().eps;
	vector_f arc;
	typename TFJVOFGen_Staggered<ST,CT>::vof_curve_type& 
			vc=this->vof.vof_curve(&arc);
	for (size_t i=0; i<vc.seg_count(); i++)
	{
		ds=vc.arcs().seg_len(i);
		// make sure don't get zero width segment
		if (ds<eps) continue;
		// detect which seg is belongs to curve
		if ((arc[i]<0) || (arc[i+1]<0)) continue;

		// setup sigma local polynomial
		sga=_sigInterp->splint(arc[i]);
		sgb=_sigInterp->splint(arc[i+1]);
		sg.set_coef(sgb-sga,sga);
		sg1=sg>>1;

		p=vc.interps().x().get_seg_poly(i);
 		q=vc.interps().y().get_seg_poly(i);
		p1=p>>1;	p2=p1>>1;
		q1=q>>1;	q2=q1>>1;

        pl=(p*p1*sg1+sg*(p2*p-q1*q1))<<1;
		ust+=pl.eval_diff(0.0,1.0)/ds;
		pl=(p*q*q1)<<1;
		ust-=Bo*pl.eval_diff(0.0,1.0);

        pl=(p*q1*sg1+sg*(p*q2+p1*q1))<<1;
		vst+=pl.eval_diff(0.0,1.0)/ds;
		pl=(q*p*p1)<<1;
		vst+=Bo*pl.eval_diff(0.0,1.0);
	}
}

template <class ST,	class CT>
void TFJVOFGen_Staggered_ST<ST,CT>::get_staggered_st(size_t m, size_t n,
						typename TFJVOFGen_Staggered<ST,CT>::staggered_type stag_loc,
						float_t& ust, float_t& vst)
{
	ust=0.0;	vst=0.0;
	int nxseg,nyseg;
	for (size_t cell=0; cell<7; cell+=2)
	{
		staggered_mn(m,n,stag_loc,cell,nxseg,nyseg);
		nxseg=bound(nxseg,0,(int)this->naxis.x().seg_count()-1);
		nyseg=bound(nyseg,0,(int)this->naxis.y().seg_count()-1);
		float_t u,v;
		if (uniform_st())
			get_cell_st(nxseg,nyseg,u,v);
		else
			get_cell_stg(nxseg,nyseg,u,v);
		ust+=u;	vst+=v;
	}
}
      
}	// end of namespace
