#include "fjapp_vof.h"
#include "fjlib_vecmat.h"

namespace fjlib {

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::insert_pair(size_t key, float_t v)
{
	// if key valid then append, otherwise insert
	seg_map_type::iterator iter=seg_map.find(key);	
	if (iter!=seg_map.end())	
		push_back((*iter).second,v);
	else
	{
		vector_f a(1); a[0]=v;
		seg_map.insert(std::make_pair(key,a));
	}
}

template <class ST,	class CT>
bool TFJVOFGen<ST,CT>::_root_unique(float_t root, vector_f& vec)
{
	if (vec.size()<1) return true;
	bool found=false;
	for (size_t i=0; i<vec.size(); i++)
		if (fabs(root-vec[i])<eps) { found=true; break; }
	return !found;
}

template <class ST,	class CT>
bool TFJVOFGen<ST,CT>::_root_unique(float_t root, int m, int n)
{
	bool found=false;
	for (size_t loc=1; loc<8; loc+=2)
	{
		size_t key=stoid(m,n,loc);
        seg_map_type::iterator iter=seg_map.find(key);	
		if (iter!=seg_map.end())	
			if (!_root_unique(root,(*iter).second)) {
				found=true; break;
			}
	}
	return !found;
}

template <class ST,	class CT>
vector_n TFJVOFGen<ST,CT>::_troot(float_t root,	scale_type& scl)
{
	vector_n minor;
	for (int u=-1; u<=1; u++)
	{
//		float_t xx=itp.splint(root+eps*u);
		float_t xx=root+eps*u;
//		if (scl.inside(xx))
		if (between(xx,scl.min()-eps,scl.max()+eps))
		{
			size_t i=scl.locate_seg(xx);
			bool ad=true;
			if (minor.size()>0)
				for (size_t v=0; v<minor.size(); v++)
					if (minor[v]==i) { ad=false; break; }
			if (ad) push_back(minor,(int)i);
		}
	}
	return minor;
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::vof_gen()
{
	seg_map.clear();
	// storage for roots
	vector_f roots;
	size_t m,n,loc;
	vector_n minor;

	size_t xsegs=xseg_count(),			// x=
		   ysegs=yseg_count();			// y=
	typename curve_type::interp_list_type& 
			interps=curve->interps();	// helper

	loc=bptS;
	// figure out y=c intersects;
	for (size_t j=0; j<ysegs+1; j++)
	{
		n=j; if (j==ysegs) { n--; loc=bptN; }
		if (interps.y().splinx(
				axis->y().value_at(j),&roots,eps))
			for (size_t k=0; k<roots.size(); k++)
			{
//				minor=_troot(roots[k],interps.x(),axis->x());
				minor=_troot(interps.x().splint(roots[k]),axis->x());
				for (size_t v=0; v<minor.size(); v++)
//					if (_root_unique(roots[k],minor[v],n))
						insert_pair(stoid(minor[v],n,loc),roots[k]);
			}
	}

	loc=bptW;
	// figure out x=c intersects;
	for (size_t i=0; i<xsegs+1; i++)
	{
		m=i; if (m==xsegs) { m--; loc=bptE; }
		if (interps.x().splinx(
				axis->x().value_at(i),&roots,eps))
			for (size_t k=0; k<roots.size(); k++)
			{
//				minor=_troot(roots[k],interps.y(),axis->y());
				minor=_troot(interps.y().splint(roots[k]),axis->y());
				for (size_t v=0; v<minor.size(); v++)
//					if (_root_unique(roots[k],m,minor[v]))
						insert_pair(stoid(m,minor[v],loc),roots[k]);
			}
	}
}

template <class ST,	class CT>
//TFJVOFGen<ST,CT>::point_type
TFJCurve::point_type TFJVOFGen<ST,CT>::box_pt(size_t i)
{
	if (vrts[i]<0)
		return ctopos(ii,jj,vloc[i]);
	else
		return curve->splint(vrts[i]);
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::sort_box()
{
	// make root unique
	size_t i=0;
	while (i<vrts.size()-1)
	{
		bool rpt=false;
		size_t j=i+1;
		while (j<vrts.size())
		{
			if (fabs(vrts[i]-vrts[j])<eps)	{ rpt=true; break; }
			j++;
		}
		if (rpt)
		{
			// delete un-unique one
			if (j<vrts.size()-1)
				for (size_t k=j; k<vrts.size()-1; k++)
				{ vrts[k]=vrts[k+1]; vloc[k]=vloc[k+1]; }
			vrts.resize(vrts.size()-1);
			vloc.resize(vloc.size()-1);
		}
		else 
			i++;
	}

	// sort by arc length
	for (size_t p=0; p<vrts.size()-1; p++)
		for (size_t q=p+1; q<vrts.size(); q++)
		{
			bool swp=false;
			if (vrts[p]>vrts[q]) swp=true;
			if ((vrts[p]==vrts[q]) && (vloc[p]<vloc[q])) swp=true;
			if (swp)
			{
				swap(vrts[p],vrts[q]);
				swap(vloc[p],vloc[q]);
			}
		}
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::tag_box()
{
	vrts.resize(0); vloc.resize(0);
	// Add roots for 4 segments of  the box
	for (size_t k=0; k<4; k++)
	{
		int kk=2*k+1;
		seg_map_type::iterator iter=seg_map.find(stoid(ii,jj,kk));
		if (iter!=seg_map.end())	// has roots on that seg
		{
			vector_f& rts=(*iter).second;
			for (size_t p=0; p<rts.size(); p++) 
			{
				push_back(vrts,rts[p]);
				push_back(vloc,kk);
			}
		}
	}
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::corner_box()
{
	vector_f urts;
	vector_n uloc;
	// Insert box corner in-between
	for (size_t i=0; i<vrts.size(); i++)
	{
		push_back(urts,vrts[i]);
		push_back(uloc,vloc[i]);

		size_t j=inc(i);
		if (vloc[i]!=vloc[j])		
		{
			int k=vloc[i];
			do
			{
				k--; if (k<0) k=7;
				if (k==vloc[j]) break;
				if (k/2*2==k)
				{
					push_back(urts,(float_t)-1.0);
					push_back(uloc,k);
				}
			} while (true);
		}
	}
	vrts=urts;
	vloc=uloc;
}

template <class ST,	class CT>
size_t TFJVOFGen<ST,CT>::inc(size_t i)
{
	size_t j=i+1; if (j==vrts.size()) j=0;
	return j;
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::insert_pt(const point_type& pt,
								 const float_t arc)
{
	for (size_t i=0; i<2; i++)
		push_back(ppos.at(i),pt.at(i));
	push_back(vof_arc,arc);
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::insert_line(size_t in, size_t jn)
{
	pos_vec_type pos(2);
	// insert points to TFJCurve_Line and TFJCurve_Poly
	int k=in;
	do
	{
		point_type p=box_pt(k);
		for (size_t i=0; i<2; i++)
			push_back(pos.at(i),p.at(i));
		if (k!=in) insert_pt(p,vrts[k]);
		k=inc(k);
	} while (k!=inc(jn));
	TFJCurve_Line cl;
	cl.set_data(&pos.at(0),&pos.at(1));

	// insert segment to TFJCurve_Poly
	for (size_t i=0; i<cl.arcs().seg_count(); i++)
	{
		for (size_t j=0; j<2; j++)
			pcurve.interp(j).polys().push_back(
						cl.interp(j).get_seg_poly(i));
	}
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::_push_poly_value(size_t xy, const TFJPolynomial &p)
{
	pcurve.interp(xy).polys().push_back(p);
/*
	if (p.degree()>1)
		pcurve.interp(xy).polys().push_back(p);
	else // hack to make TFJVOFGen_ST calculation work for linear curve
	{
		TFJPolynomial q=p;
		q.set_degree(2);
		pcurve.interp(xy).polys().push_back(q);
	}
*/
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::insert_poly(size_t in, size_t jn)
{
	// both i, and j have to be non=cornered 
	size_t iseg=curve->arcs().locate_seg(vrts[in]);
	size_t jseg=curve->arcs().locate_seg(vrts[jn]);
	float_t irel=curve->arcs().seg_relative(iseg,vrts[in]);
	float_t jrel=curve->arcs().seg_relative(jseg,vrts[jn]);

	size_t un;
	if (iseg<jseg)
	{
		un=upper_pt(iseg);
		insert_pt(curve->pt_at(un),
					curve->arcs().value_at(un));
		// add partial curve 
		for (size_t i=0; i<2; i++)
		{
			TFJPolynomial p=curve->interp(i).get_seg_poly(iseg);
			p%=irel;	p&=1-irel;
//			pcurve.interp(i).polys().push_back(p);
			_push_poly_value(i,p);
		}
		// add middle full pieces
		if (iseg<jseg-1)
		{
			for (size_t k=iseg+1; k<jseg; k++)
			{
				un=upper_pt(k);
				insert_pt(curve->pt_at(upper_pt(k)),
							curve->arcs().value_at(un));
				for (size_t i=0; i<2; i++)
//					pcurve.interp(i).polys().push_back(
//						curve->interp(i).get_seg_poly(k));
				_push_poly_value(i,
						curve->interp(i).get_seg_poly(k));
			}
		}
		insert_pt(box_pt(jn),vrts[jn]);
		// add partial curve 
		for (size_t i=0; i<2; i++)
		{
			TFJPolynomial p=curve->interp(i).get_seg_poly(jseg);
			p&=jrel;
//			pcurve.interp(i).polys().push_back(p);
			_push_poly_value(i,p);
		}
	}	
	else		// iseg==jseg
	{
		// add small piece
		insert_pt(box_pt(jn),vrts[jn]);
		// add partial curve 
		for (size_t i=0; i<2; i++)
		{
			TFJPolynomial p=curve->interp(i).get_seg_poly(iseg);
			p%=irel;	p&=jrel-irel;	
//			pcurve.interp(i).polys().push_back(p);
			_push_poly_value(i,p);
		}
	}
}

template <class ST,	class CT>
bool TFJVOFGen<ST,CT>::inside_box(size_t in, size_t jn)
{
	float_t hs=(vrts[in]+vrts[jn])/2.0;
	point_type p=curve->splint(hs);
	return (between(p.x(),axis->x().seg_lower_value(ii)-eps,
						axis->x().seg_upper_value(ii)+eps) &&
			between(p.y(),axis->y().seg_lower_value(jj)-eps,
						axis->y().seg_upper_value(jj)+eps));
/*
	return ((axis->x().inside_seg(p.x(),ii)) &&
			(axis->y().inside_seg(p.y(),jj)));
*/
}

template <class ST,	class CT>
void TFJVOFGen<ST,CT>::link_box()
{
	// clear the node pos list
	ppos.redim(2);
	pcurve.interps().redim(2);
	for (size_t i=0; i<2; i++)
	{
		ppos.at(i).resize(0);
		pcurve.interp(i).polys().resize(0); 
	}

	// clear vof_arc
	vof_arc.resize(0);
	
	point_type n1,n2;
	size_t ps=0;
	// find the first non-corner point
//	for (size_t i=0; i<vrts.size(); i++)
//		if (vrts[i]>0) { ps=i; break; }
	insert_pt(box_pt(ps),vrts[ps]);
	
	int p=ps,q;
	do 
	{
		// find the next non-corner point
		q=inc(p);
		while (vrts[q]<0)
		{ q=inc(q); }
		if (q==ps) break;
		// check if this piece inside the box
		if (inside_box(p,q))
			insert_poly(p,q);
		else
			insert_line(p,q);
		p=q;
	} while (true);
	insert_line(p,q);

	pcurve.set_data(&ppos.at(0),&ppos.at(1));
}

template <class ST,	class CT>
bool TFJVOFGen<ST,CT>::is_colored(size_t in, size_t jn)
{
	ii=in; jj=jn;
	tag_box();
	if (vrts.size()>1)
		sort_box();
	return (vrts.size()>1);
}

template <class ST,	class CT>
//TFJVOFGen<ST,CT>::vof_curve_type& 
TFJCurve_Poly& TFJVOFGen<ST,CT>::vof_curve(vector_f* arc_info)
{
	// Check is_colored() before calling vof_curve() 
	corner_box();
	link_box();
	if (arc_info!=NULL) *arc_info=vof_arc;
	return pcurve;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
template <class ST,	class CT>
void TFJVOFGen_Full<ST,CT>::get_ycenter_roots_xseg(size_t ys, vector_n& xrts)
{
	size_t xsegs=xseg_count(),			// x=
		   ysegs=yseg_count();			// y=
	typename curve_type::interp_list_type&
		interps=curve->interps();		// helper

	xrts.resize(0);
	// get roots
	vector_f rts;
	if (!interps.y().splinx(
			axis->y().seg_center_value(ys),&rts,eps))
		return;

	// make it unique, different one here
	if (rts.size()>1)
	{
		vector_f tmp(1);
		tmp[0]=rts[0];
		if (rts.size()>2)
			for (size_t i=1; i<rts.size(); i++)
			{
				bool add=true;
				for (size_t j=0; j<i; j++)
					if (fabs(rts[j]-rts[i])<eps) 
					{ add=false; break; }
				if (add) push_back(tmp,rts[i]);
			}
		rts=tmp;
	}

	// convert roots to xpos 
	for (size_t i=0; i<rts.size(); i++)
		rts[i]=interps.x().splint(rts[i]);

	// sort by xpos
	if (rts.size()>1)
		for (size_t p=0; p<rts.size()-1; p++)
			for (size_t q=p+1; q<rts.size(); q++)
				if (rts[p]>rts[q]) 
					swap(rts[p],rts[q]);

	// convert xpos to xseg
	xrts.resize(rts.size());
	for (size_t i=0; i<rts.size(); i++)
		xrts[i]=axis->x().locate_seg(rts[i]);
}

template <class ST,	class CT>
float_t TFJVOFGen_Full<ST,CT>::get_cell_vof(size_t m,size_t n)
{
	if (is_colored(m,n))
	{
		float_t sum=0.0;
		vof_curve_type& vc=vof_curve();
		for (size_t i=0; i<vc.arcs().seg_count(); i++)
		{
			TFJPolynomial p=vc.interps().x().get_seg_poly(i),
						  q=vc.interps().y().get_seg_poly(i);
			if (!cylinder)
				p=(p*(q>>1))<<1;
			else
				p=(p*p*(q>>1))<<1;
			float_t sm=(p|1.0)-(p|0.0);
			sum+=sm;
		}
		if (!cylinder)
			return  -sum/get_cell_area(m,n);
		else
			return -sum/2.0/axis->x().seg_center_value(m)/
							get_cell_area(m,n);
/*		
		if (!cylinder)
			return  -sum/axis->x().delta()
						/axis->y().delta();
		else
			return -sum/2.0/axis->x().seg_center_value(m)
						/axis->x().delta()
						/axis->y().delta();
*/
	}
	else return void_value;
}

template <class ST,	class CT>
void TFJVOFGen_Full<ST,CT>::calc_vofs(matrix_f &mat)
{
	size_t m=xseg_count(),
			n=yseg_count();
	mat.resize(m,n);

	for (size_t i=0; i<m; i++)
		for (size_t j=0; j<n; j++)
			mat(i,j)=get_cell_vof(i,j);
}

template <class ST,	class CT>
void TFJVOFGen_Full<ST,CT>::fill_vofs(matrix_f &mat)
{
	size_t xn=xseg_count();
	size_t yn=yseg_count();

	size_t cn=curve->pt_count();
	size_t ymax=axis->y().locate_seg(curve->y(0));
	size_t ymin=axis->y().locate_seg(curve->y(cn-1));

	size_t j=0;
	while (j<yn)
	{
		vector_n xrts;
		get_ycenter_roots_xseg(j,xrts);
		if (xrts.size()<1) { j++; continue; }

		bool _bound=((j<=ymax) && (j>=ymin));
		bool _val=_bound;
		for (size_t i=0; i<xrts.size(); i++)
		{
			size_t m,n;
			// handle the first node
			if (i==0)
			{
				if (_bound) { m=0; n=xrts[0]; }
				else continue;
			}
			else { m=xrts[i-1]; n=xrts[i]; }
			if (_val)
				for (size_t k=m; k<=n; k++)
					if (mat(k,j)==void_value) mat(k,j)=1.0;
			_val=!_val;
		}
		j++;
	}

	// replace void_value with 0
	for (size_t m=0; m<mat.size1(); m++)
		for (size_t n=0; n<mat.size2(); n++)
			if (mat(m,n)==void_value) mat(m,n)=0.0;
}
   
}	// end of namespace