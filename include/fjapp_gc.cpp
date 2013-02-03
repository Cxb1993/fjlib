#include "fjapp_gc.h"
#include "fjlib_vecmat_print.h"
#include "fjlib_cio.h"
#include "fjlib_blas.h"
#include "fjlib_polyroots.h"

namespace fjlib {

template <class ST,class CT>
void TFJGCGen<ST,CT>::gen_unknown_map()
{
	_bnodes.resize(0);
	_ghosts.resize(0);
	size_t xn=_umat.size1(),
			yn=_umat.size2();
			
	vector_n xx,yy;
	set_vector(xx,0,0,-1,1);
	set_vector(yy,-1,1,0,0);
	npos_type a(2);
	for (size_t j=0; j<yn; j++)
		for (size_t i=0; i<xn; i++)
		{
			if (_umat(i,j)!=0) continue;
			// find a vacant site
			bool is_ghost=false;
			for (size_t k=0; k<4; k++)
			{
				int m=xx[k]+i, n=yy[k]+j;
				if (!between(m,0,(int)xn)) continue;
				if (!between(n,0,(int)yn)) continue;
				// find a nearby non-vacant site
				if (_umat(m,n)>0) is_ghost=true;
				if (_umat(m,n)==1) {
					a.x()=m; a.y()=n; _bnodes.push_back(a);
					_umat(m,n)=2;
				}
			}
			if (is_ghost) {
				a.x()=i; a.y()=j; _ghosts.push_back(a);
				_umat(i,j)=-1;
			}
		}
}

template <class ST,class CT>
float_t TFJGCGen<ST,CT>::_dist(float_t x1, float_t y1,
								float_t x2, float_t y2)
{
	float_t dx=x1-x2, dy=y1-y2;
	return sqrt(dx*dx+dy*dy);
}
								
template <class ST,class CT> inline
size_t TFJGCGen<ST,CT>::_bound_max(size_t v, size_t m)
{
	if (m==0) throw "max cant be zero in fjapp_gc.cpp";
	if (v<m) return v;
	else return (size_t)(m-1);
}

template <class ST,class CT> inline
size_t TFJGCGen<ST,CT>::_bound_min(size_t v, size_t dv)
{
	if (v>=dv) return v-dv;
	else return 0;
}

template <class ST, class CT>
bool TFJGCGen<ST,CT>::_get_range(float_t x, float_t y, int n,
						size_t& mx1, size_t& mx2, 
						size_t& my1, size_t& my2)
{
	// locate x,y
	size_t nx,ny;
	nx=_axis->x().locate_seg(x);
	ny=_axis->y().locate_seg(y);
	// estimate range
	int dn=(int)(sqrt((float)n)+0.5);
	
	if (dn==0) dn=1;
	// check range
	mx1=_bound_min(nx,dn);
	mx2=_bound_max(nx+dn,_axis->x().seg_count());
	if (mx2<mx1) return false;
	my1=_bound_min(ny,dn);
	my2=_bound_max(ny+dn,_axis->y().seg_count());
	if (my2<my1) return false;
	
	return true;
}

template <class ST, class CT>
bool TFJGCGen<ST,CT>::_sort_vec(const vector_f& v, size_t n,
							vector_n& od)
{
	if (n>v.size()) throw "can't sort too many elements in fjapp_gc.cpp";
	vector_f vc=v;
	od.resize(v.size());
	for (size_t i=0; i<od.size(); i++)
		od[i]=i;
	// sort
	size_t kn=n;
	if (kn==0) kn=v.size();
	float_t _min,tmp;
	for (size_t i=0; i<kn; i++)
	{
		if (i>=vc.size()-1) continue;	// can't do the last one
		_min=1e30;
		size_t k=i;
		for (size_t j=i; j<vc.size(); j++)
			if (vc[j]<_min) { _min=vc[j]; k=j; }
//		print_var("k",k);
		if (k!=i) { 
			tmp=vc[i]; vc[i]=vc[k]; vc[k]=tmp;
			tmp=od[i]; od[i]=od[k]; od[k]=tmp;
		}
	}
	return true;
}

								
template <class ST,class CT>
void TFJGCGen<ST,CT>::search_grid_nodes(float_t x, float_t y,
								node_vec_type& nvec, size_t n)
{
	if (n<1) return;
	// make a list of grid nodes
	size_t mx1,mx2,my1,my2;
	// get range
	if (!_get_range(x,y,n,mx1,mx2,my1,my2))
		throw "cant range enough grid nodes in fjapp_gc.cpp";
//	cout << "(" << x << "," << y << ") "
//		<< mx1 << "-" << mx2 << "  " << my1 << "-" << my2 << endl;
	size_t gnx, gny,gn;
	// make a list of grid nodes
	vector_f dist;
	vector_n vx,vy; 
	while (true) 
	{
		gnx=mx2-mx1+1;	gny=my2-my1+1;	gn=gnx*gny;
		if (gn<n) throw "_get_range() doesnt return enough nodes in fjapp_gc.cpp";
//		dist.resize(gn);
		vx.resize(0); vy.resize(0);
		dist.resize(0);
//		size_t k=0,ku=0;
		float xg,yg;
		for (size_t p=(size_t)mx1; p<=(size_t)mx2; p++)
		{
			xg=_axis->x().seg_center_value(p);
			for (size_t q=(size_t)my1; q<=(size_t)my2; q++)
			{
				yg=_axis->y().seg_center_value(q);
				if (_umat(p,q)>0) 
				{
					push_back(vx,(int)p);
					push_back(vy,(int)q);
					push_back(dist,_dist(xg,yg,x,y));
					if (p==0) // take care symemtric here
					{
						push_back(vx,-1);
						push_back(vy,(int)q);
						push_back(dist,_dist(-xg,yg,x,y));
					}
				}
/*	
				{ dist[k]=_dist(xg,yg,x,y); ku++; }
				else dist[k]=1e30;		// exclude the known
				k++;
*/				
			}
		}
//		if (ku>=n) break;	// found
		if (dist.size()>=n) break; // found
		mx1=_bound_min(mx1,1);	mx2++;
		if ((mx2-mx1+1)>_axis->x().seg_count()) 
			throw "cant find enough nodes on whole domain x in fjapp_gc.cpp";
		mx2=_bound_max(mx2,_axis->x().seg_count());
		my1=_bound_min(my1,1); my2++;
		if ((my2-my1+1)>_axis->y().seg_count()) 
			throw "cant find enough nodes on whole domain y in fjapp_gc.cpp";
		my2=_bound_max(my2,_axis->y().seg_count());
	}
//	cout << "mx: " << mx1 << " " << mx2 << endl;
//	cout << "my: " << my1 << " " << my2 << endl;
//	print_var("dist",dist);
	// sort the list and pick the first nth
	vector_n od;
	if (!_sort_vec(dist,n,od))
		throw "cant sort enough grid nodes in fjapp_gc.cpp";
	// record
	size_t ind;
	int ii,jj;
	for (size_t k=0; k<n; k++)
	{
		ind=od[k];
		ii=vx[ind];	jj=vy[ind];
//		ii=ind/gny; jj=ind-ii*gny;
//		ii+=mx1;	jj+=my1;
		node_type nd;
		if (ii<0)		// symetry handling 
			nd.x=-_axis->x().seg_center_value(0);
		else
			nd.x=_axis->x().seg_center_value(ii);
		nd.y=_axis->y().seg_center_value(jj);
		nd.type=2;	nd.i=ii; nd.j=jj;
		nd.dist=dist[ind];
		nvec.push_back(nd);
	}
}

template <class ST,class CT>
void TFJGCGen<ST,CT>::search_surf_nodes(float_t x, float_t y,
								node_vec_type& nvec, size_t n)
{
	if (n<1) return;
	// get a list of distantce 
	size_t cn=_curve->pt_count();
	vector_f dist(cn);
	float xc,yc;
	for (size_t i=0; i<cn; i++)
	{
		xc=_curve->x(i); yc=_curve->y(i);
		dist[i]=_dist(xc,yc,x,y);
	}
	// sort this list and pick the first nth 
	vector_n od;
	if (!_sort_vec(dist,n,od))
		throw "cant sort enough curve nodes in fjapp_gc.cpp";
	// record
	size_t j;
	for (size_t k=0; k<n; k++)
	{
		node_type nd;
		j=od[k];
		nd.x=_curve->x(j); nd.y=_curve->y(j);
		nd.type=1;	nd.i=j;
		nd.dist=dist[j];
		nvec.push_back(nd);
	}
}

/*
template <class ST,class CT>
void TFJGCGen<ST,CT>::convert_npos_to_pos(const npos_vec_type& vp)
{
	size_t gn=_ghosts.size();
	vp.resize(gn);
	for (size_t i=0; i<gn; i++)
	{
		vp[i].redim(2);
		vp[i].x()=_axis->x().seg_center_value(_ghosts[i].x());
		vp[i].y()=_axis->y().seg_center_value(_ghosts[i].y());
	}
}
*/

template <class ST,class CT>
void TFJGCGen<ST,CT>::gen_env_info(const pos_vec_type& vp,
									size_t bn,size_t sn,bool sort)
{
	size_t tn=bn+sn;
	if ((tn!=3) && (tn!=6)) 
		throw "only take 3 points or 6 points interpolation in fjapp_gc.cpp";
	size_t gn=vp.size();
	_env_info->resize(gn);
	// for each ghost
	for (size_t i=0; i<gn; i++)
	{	
		node_vec_type &nv=(*_env_info)[i];
		nv.resize(0);
		search_grid_nodes(vp[i].x(),vp[i].y(),nv,bn);
		search_surf_nodes(vp[i].x(),vp[i].y(),nv,sn);
/*
		cout << "#" << i << " (" << _ghosts[i].x() 
			<< "," << _ghosts[i].y() << ")" << endl;
		for (size_t j=0; j<nv.size(); j++)
			cout << nv[j].x << "," << nv[j].y 
				<< "[" << nv[j].i << "," << nv[j].j << "]D"
				<< nv[j].dist << endl;
*/
	}
	if (!sort) return;
/* wait later to implement	
	if (vp.size()<2) return;
	// sort by dist
	for (size_t i=0; i<vp.size()-1; i++)
	{
		// first the minimum
		float_t dist=vp[i].dist;
		size_t k=i;
		for (size_t j=i+1; j<vp.size(); j++)
			if (vp[j].dist<dist) { k=j; dist=vp[j].dist; }
		if (k!=i)
		{
			node_type tmp=vp[k];
			vp[k]=vp[i]; vp[i]=tmp;
		}
	}
*/
}

/*
template <class ST,class CT>
void TFJGCGen<ST,CT>::gen_env_info(size_t type, size_t bn, size_t sn)
{
	if (type==0) 	// not predefined, call the other one
		throw "please call orginal version of gen_env_info in fjapp_gc.cpp";
	pos_vec_type vp;
	if (type==1)	// on surf node
	{
		vp.resize(_curve->pt_count());
		for (size_t i=0; i<vp.size(); i++)
		{
			vp[i].redim(2);
			vp[i].x()=_curve->x(i);	vp[i].y()=_curve->y(i);
		}
	}
	if (type==2)	// on ghost grid node
	{
		vp.resize(_ghosts.size());
		for (size_t i=0; i<vp.size(); i++)
		{
			vp[i].redim(2);
			vp[i].x()=_axis->x().seg_center_value(_ghosts[i].x());
			vp[i].y()=_axis->y().seg_center_value(_ghosts[i].y());
		}
	}
	gen_env_info(vp,bn,sn);
}
*/

template <class ST,class CT>
void TFJGCGen<ST,CT>::fill_env_nodes(const matrix_f& vb,
									const vector_f& vs)
{
	for (size_t i=0; i<_env_info->size(); i++)
	{
		node_vec_type& nv=(*_env_info)[i];
		for (size_t j=0; j<nv.size(); j++)
		{
			if (nv[j].type==1)	// on curve
				nv[j].v=vs[nv[j].i];
			if (nv[j].type==2)	// on bulk
				if (nv[j].i<0)
					nv[j].v=vb(0,nv[j].j);	// symmetry
				else
					nv[j].v=vb(nv[j].i,nv[j].j);
		}
	}
}

template <class ST,class CT>
void TFJGCGen<ST,CT>::gen_interp_coef(const vector_n& src,
									matrix_f& coef)
{
	size_t nn=src.size();
	coef.resize(_env_info->size(),nn);
	if (nn==6) {		// 6 point interpolation
	// still got chance to fail if two nodes has identical locations, rare
// matrix inverse	
	matrix_f lc(nn,nn);
	vector_f lv(nn),la(nn);
	// for each ghost node
	for (size_t i=0; i<_env_info->size(); i++)
	{
		node_vec_type& nv=(*_env_info)[i];
		// assemble local matrix 
		size_t kk;
		float_t x,y;
		for (size_t j=0; j<nn; j++)
		{
			lc(j,0)=1.0;
			kk=src[j];	x=nv[kk].x;	y=nv[kk].y;
			lc(j,1)=x;		lc(j,2)=y;
			lc(j,3)=sqr(x);	lc(j,4)=x*y;
			lc(j,5)=sqr(y);
			lv[j]=nv[kk].v;
		}
		// solve coef lc*la=lv
		solve_axb(lc,lv,la);
	/*
		print_vec("lv",lv);
		print_mat("lc",lc);
		print_vec("coef_inv",la);
	*/
		for (size_t j=0; j<nn; j++)
			coef(i,j)=la[j];
	}
	} else {			// 3 point interpolation
	// 3 point matrix inverse has sigular problem when coordinates same in one direction
/* analytic solution */
	vector_f x(nn),y(nn),t(nn),a(nn);
	matrix_f m(nn,nn);
	float_t a2;
	size_t i=0,j=1,k=2;		
	for (size_t p=0; p<_env_info->size(); p++)
	{
		node_vec_type& nv=(*_env_info)[p];
		size_t kk;
		for (size_t q=0; q<nn; q++)
		{
			kk=src[q];
			x[q]=nv[kk].x; y[q]=nv[kk].y;	t[q]=nv[kk].v;
//			cout << x[q] << "," << y[q] << ":" << t[q] << endl;
		}
		a2=x[i]*(y[j]-y[k])+x[j]*(y[k]-y[i])+x[k]*(y[i]-y[j]);
		if (fabs(a2)<1e-12) {
			// it's a hack here for my purpose, should throw in general
			reset_vector(a,0.0);	a[0]=t[i];
		}
		else {
			m(0,0)=x[j]*y[k]-x[k]*y[j];
			m(0,1)=x[k]*y[i]-x[i]*y[k];
			m(0,2)=x[i]*y[j]-x[j]*y[i];
			m(1,0)=y[j]-y[k];
			m(1,1)=y[k]-y[i];
			m(1,2)=y[i]-y[j];
			m(2,0)=x[k]-x[j];
			m(2,1)=x[i]-x[k];
			m(2,2)=x[j]-x[i];
//			print_mat("m",m);
//			print_var("t",t);
//			print_var("a2",a2);
			a=(1.0/a2)*ublas::prod(m,t);
		}
//		print_vec("coef_an",a);
		for (size_t q=0; q<nn; q++)
			coef(p,q)=a[q];
	}
	}
}

template <class ST,class CT>
float_t TFJGCGen<ST,CT>::interp_xy(size_t index,
									float_t x, float_t y,
									const matrix_f& coef)
{
	size_t nn=coef.size2();
	vector_f c(nn);
	for (size_t i=0; i<nn; i++)
		c[i]=coef(index,i);
	float_t v=c[0]+c[1]*x+c[2]*y;
	if (nn==6)
		v+=c[3]*sqr(x)+c[4]*x*y+c[5]*sqr(y);
	return v;
}

template <class ST,class CT>
void TFJGCGen<ST,CT>::gc_gen()
{
	gen_unknown_map(); 
	_prepare_pos();
	gen_mirrors();
}

template <class ST,class CT>
void TFJGCGen<ST,CT>::_prepare_pos()
{
	size_t i;
	// ghost pos
	_gpos.resize(_ghosts.size());
	i=0;
	while (i<_gpos.size()) {
		_gpos[i].redim(2);
		_gpos[i].x()=_axis->x().seg_center_value(_ghosts[i].x());
		_gpos[i].y()=_axis->y().seg_center_value(_ghosts[i].y());
		i++;
	}
	// curve pos
	_spos.resize(_curve->pt_count());
	for (i=0; i<_spos.size(); i++)
	{
		_spos[i].redim(2);
		_spos[i].x()=_curve->x(i);
		_spos[i].y()=_curve->y(i);
	}
}

/*
template <class ST,class CT>
void TFJGCGen<ST,CT>::gen_mirror_ghosts()
{
	size_t gn=_ghosts.size();
	_mgpos.resize(gn);
	_ms_npos.resize(gn);
	node_vec_type gv;
	size_t j;
	for (size_t i=0; i<gn; i++)
	{
		// search for the closest node on surface
		gv.resize(0);
		search_surf_nodes(_gpos[i].x(),_gpos[i].y(),gv,1);
		j=_ms_npos[i]=gv[0].i;
		// make mirror
		_mgpos[i].redim(2);
		_mgpos[i].x()=2.0*_spos[j].x()-_gpos[i].x();
		_mgpos[i].y()=2.0*_spos[j].y()-_gpos[i].y();
		cout << _mgpos[i].x() << " " << _mgpos[i].y() << endl;
	}
}
*/

template <class ST,class CT>
void TFJGCGen<ST,CT>::gen_mirrors()
{
	size_t gn=_ghosts.size();
	_mpos.resize(gn);
	_marc.resize(gn);
	if (gn==0) return;
//	_mgpos.resize(gn);
	
	// create a lcurve due to polyroots can't solve higher than order of 6
	TFJCurve_Line cv;
	typedef typename curve_type::point_list_type pt_type;
	const pt_type &t=_curve->get_points();
	cv.set_data(t.x(),t.y());
	cv.calc_arcs();
	size_t nseg=cv.seg_count();
	std::vector<TFJPolynomial> px(nseg),py(nseg);
	for (size_t i=0; i<nseg; i++)
	{
		px[i]=cv.interps().x().get_seg_poly(i);
		py[i]=cv.interps().y().get_seg_poly(i);
	}

	// for each ghost node, find the mirror node on the surface
	TFJPolyRoots pr;
	vector_f rts;
	TFJPolynomial pxx,pyy,p;
	size_t rn;
	for (size_t i=0; i<gn; i++)
	{
		node_vec_type nv;
		for (size_t j=0; j<nseg; j++)
		{
			pxx=px[j];	pxx[0]-=_gpos[i].x();
			pyy=py[j];	pyy[0]-=_gpos[i].y();
			p=(px[j]>>1)*pxx+(py[j]>>1)*pyy;
//			cout << p << endl;
			size_t n=pr.solve(p);
//			cout << i << " " << j << " " << rts.size() << endl;
			if (n<1) continue;
			rts.resize(0);
			pr.real_roots(rts);
//			cout << i << "," << j << ")" << rts << endl;
			if (rts.size()<1) continue;
			else {
				size_t k=0;
				while (k<rts.size())
				{
					if (between(rts[k],-0.1,1.1)) break;
					k++;
				}
				if (k<rts.size())
				{
					node_type a;
					a.i=j;			// record segment id
					float_t sc=bound(rts[k],0.0,1.0);
					float_t ss=_curve->arcs().seg_balance(j,sc);
					typename curve_type::point_type b=_curve->splint(ss);
					a.x=b.x();	a.y=b.y();
					a.dist=ss;
					a.type=1;		// curve sega
					a.v=rts[k];
					nv.push_back(a);
				}
			}
		}
		_mpos[i].redim(2);
		// pick the nearest one
		if (nv.size()<1)
		{
			cout << "warning: can't find normal node for ghost #" 
				<< i << "(" << _gpos[i].x() << "," << _gpos[i].y()
				<< endl;
//			throw "can't find mirror in fjapp_gc.cpp";
			_mpos[i].x()=_gpos[i].x();
			_mpos[i].y()=_gpos[i].y();
			_marc[i]=-1;
		} else {
			float_t dist=1e30,tmp;
			size_t kk=0;
			for (size_t k=0; k<nv.size(); k++)
			{
				tmp=sqrt(sqr(nv[k].x-_gpos[i].x())+
						sqr(nv[k].y-_gpos[i].y()));
				if (tmp<dist) { dist=tmp; kk=k; }
			}
			_mpos[i].x()=nv[kk].x;	_mpos[i].y()=nv[kk].y;
			_marc[i]=nv[kk].dist;
/*
			if (!between(nv[kk].v,0.0,1.0))
				cout << "warning: ghost node " << i 
					<< " intercept with roots "
					<< nv[kk].v << endl;
*/
		}
	}
}

}	// end of namespace

