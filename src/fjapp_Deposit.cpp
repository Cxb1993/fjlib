#include "fjapp_SphereAdsorp.h"
#include "fjlib_vecmat_print.h"

namespace fjlib {

#ifdef _Petsc_KSP
void TFJSphereAdsorp::initialize_ksp()
{
	_cksp.set_data(&_cap,&_caw,&_cae,
					&_cas,&_can,&_csu,&_conc);
	_cksp.set_mask(&_gc.get_umat());
	_cksp.set_tolerance(_params.conc_acc);
}
#endif


void TFJSphereAdsorp::initialize()
{
	params_type &p=_params;
	
	// create a geometry coordinates, half spherical x>0
	float_t dlen=1.0/p.npl;
	int nn=int(p.size*p.npl+0.5);	dlen=p.size/nn;
	_axis.redim(2);
	_axis.x().set_range(0.0,dlen,nn);
	_axis.y().set_range(-p.size,dlen,nn*2);
	
	// intialize bulk property //	
	_conc.resize(xn(),yn());

	// create a half spherical surface
	size_t nc=M_PI*p.radius*p.npl+1;
	if (nc<5) throw "come on, give more nodes in fjapp_SphereAdsorp.cpp";
	_x.resize(nc); 		_y.resize(nc);
	_xsp.resize(nc);	_ysp.resize(nc);
	_cn1.resize(nc);
	float_t dsita=M_PI/(nc-1);
	float_t ac;
	for (size_t i=0; i<nc; i++)
	{
		ac=dsita*i;
		_x[i]=p.radius*std::sin(ac);
		_y[i]=p.radius*std::cos(ac);
		_xsp[i]=_y[i];	
		_ysp[i]=-_x[i];
		_cn1[i].redim(2);
		_cn1[i].x()=_x[i]-_ysp[i]/p.npl;
		_cn1[i].y()=_y[i]+_xsp[i]/p.npl;
	}
	_x[0]=-1e-12;
	_x[nc-1]=-1e-12;

	if (_params.debug) {
	cout << "_x,_y:" << endl;
	for (size_t i=0; i<_x.size(); i++)
		cout << i << " " << _x[i] << "," << _y[i] << endl;
	cout << "_xn1,_yn1:" << endl;
	for (size_t i=0; i<_x.size(); i++)
		cout << i << " " << _cn1[i].x() << "," << _cn1[i].y() << endl;
	}
//
	// initialize surf
	_surf.set_data(&_x,&_y);
	typedef interp_type at;
	_surf.interps().x().set_bc(at::csbtClamped,
							at::csbtClamped,1,-1);
	_surf.interps().y().set_bc(at::csbtClamped,
							at::csbtClamped,0,0);
	_surf.calc_arcs();
	_surf.spline();

	// initialize surf property
	_gama.resize(sn());
	_csub.resize(sn());
	_cfluxsub.resize(sn());

	// setup interpolation gama for mirror
	_sc_interp.set_data(_surf.arcs().get_vector(),&_csub);
	_sc_interp.set_bc(at::csbtClamped,at::csbtClamped,0,0);

	// setup vof
	_vof.set_data(&_axis,&_surf);
//	_vof.set_coordinates(true);
	_vof.get_vof().set_coordinates(true);
	_vof.vof_gen();

/* print roots info
	std::vector<vector_f>& rts=_vof.get_ycrtsx();
	for (size_t i=0; i<rts.size(); i++)
	{	
		size_t j=0;
		cout << i << ":";
		while (j<rts[i].size()) {
			cout << (rts[i])[j] << " ";
			j++;
		}
		cout << endl;
	}		
*/	
	_vof.fill_unknowns(_pre_umat);
//	print_mat("_pre_umat",_pre_umat,vmpfText | vmpfColMajor | vmpfColReverse);
//	save_mat("premat.txt",_pre_umat);
	
	// setup gc
	_gc.set_data(&_axis,&_surf,_pre_umat);
	_gc.gc_gen();

	// setup ghost mirror node
	_mgpos.resize(gn());
	_mxsp.resize(gn()); _mysp.resize(gn());
	_mn1.resize(gn());
	vector_f tmp;
	for (size_t i=0; i<gn(); i++)
	{
		gc_type::pos_type a=_gc.mirrors_pos()[i];
		gc_type::pos_type b=_gc.ghosts_pos()[i];
//
		_mgpos[i].redim(2);
		_mgpos[i].x()=2.0*a.x()-b.x();
		_mgpos[i].y()=2.0*a.y()-b.y();
//
		_surf.interps().x().splint(_gc.mirrors_arc()[i],&tmp);
		_mxsp[i]=tmp[1];
		_surf.interps().y().splint(_gc.mirrors_arc()[i],&tmp);
		_mysp[i]=tmp[1];
		
		_mn1[i].redim(2);
		_mn1[i].x()=a.x()-_mysp[i]/p.npl;
		_mn1[i].y()=a.y()+_mxsp[i]/p.npl;
	}
// 
	if (_params.debug) {
	cout << "mirror nodes:" << endl;
	for (size_t i=0; i<gn(); i++)
	{
		gc_type::pos_type a=_gc.mirrors_pos()[i];
		gc_type::pos_type b=_gc.ghosts_pos()[i];
		gc_type::pos_type c=_mn1[i];
		cout << i << " g:" << b.x() << "," << b.y() << 
					" m:" << a.x() << "," << a.y() <<
					" mn1:" << c.x() << "," << c.y()
					<< endl;
	}
	}
//
	// setup gama solver
	_gama_solver.set_data(&_gcoefB,&_gcoefA,&_gcoefC,&_gcoefR,&_gama);

	// initialize interp scheme
	initialize_interp();
//	
	if (_params.debug) {
	cout << "_senv:" << endl;
	for (size_t i=0; i<_senv.size(); i++)
	{
		cout << " " << i << ":";
		gc_type::node_vec_type& nv=_senv[i];
		for (size_t j=0; j<nv.size(); j++)
			cout << "("<< nv[j].i << "," << nv[j].j << ") ";
		cout << endl;
		for (size_t j=0; j<nv.size(); j++)
			cout << "("<< nv[j].x << "," << nv[j].y << ") ";
		cout << endl;
	}
	}
//
	// initial parallel
#ifdef _Petsc_KSP
	initialize_ksp();
#endif
}

void TFJSphereAdsorp::initialize_interp()
{
	// make surf env
	_gc.set_env(&_senv);
	_gc.gen_env_info(_cn1,3,0);		// for flux
	set_vector(_s_src,0,1,2);
	// make mirror env
	_gc.set_env(&_genv);
//	_gc.gen_env_info(_mgpos,2,1);
	_gc.gen_env_info(_mn1,3,0);		// for cg
	set_vector(_g_src,0,1,2);
}

/*
void TFJSphereAdsorp::update_fake_bc(const matrix_f& m)
{
	// neumann boundary condition on the edges
	// to faciliate calculating gradient
	int i,j;
	// up and bottom rows
	for (i=0; i<(int)xn(); i++)
	{
		j=0; 		_m(i,j)=_m(i,j+1);
		j=yn()-1; 	_m(i,j)=_m(i,j-1);
	}
	// left and right columns
	for (j=0; j<(int)yn(); j++)
	{
		i=0;		_m(i,j)=_m(i+1,j);
		i=xn()-1;	_m(i,j)=_m(i-1,j);
	}
}
*/

/*
size_t TFJSphereAdsorp::_index(size_t i, int di, bool y)
{
	int index=(int)i+di;
	if (y) return bound(index,0,(int)yn()-1);
	else return bound(index,0,(int)xn()-1);
}

void TFJSphereAdsorp::_set_grad(size_t i, size_t j, bool check)
{
	if (_gc.get_umat()(i,j)<=0) return;	// unknowns only
	if (!check)
	{
		_cgrad_x(i,j)=(_conc(i+1,j)-_conc(i-1,j))/dx()/2;
		_cgrad_y(i,j)=(_conc(i,j+1)-_conc(i,j-1))/dy()/2;
		return;
	}
	size_t ip=_index(i,1,false),
		   iq=_index(i,-1,false),
		   jp=_index(j,1,true),
		   jq=_index(j,-1,true);
//	cout << "ip " << ip << " iq " << iq << endl;
	_cgrad_x(i,j)=(_conc(ip,j)-_conc(iq,j))/dx()/2;
	_cgrad_y(i,j)=(_conc(i,jp)-_conc(i,jq))/dy()/2;
}

void TFJSphereAdsorp::update_grad_conc(bool tag2_only)
{
	size_t p=xn(),q=yn();
	_cgrad_x.resize(p,q);	_cgrad_y.resize(p,q);
	if (!tag2_only) 
	{
		size_t i,j;
		// excluding edge
		for (i=1; i<(int)xn()-1; i++)
			for (j=1; j<(int)yn()-1; j++)
					_set_grad(i,j,false);
		// top & bottom
		for (i=0; i<(int)xn(); i++)
		{
			j=0;	 	_set_grad(i,j,true);
			j=yn()-1;	_set_grad(i,j,true);
		}	
		// left & right
		for (j=1; j<(int)yn()-1; j++)
		{
			i=0;		_set_grad(i,j,true);
			i=xn()-1;	_set_grad(i,j,true);
		}
	} else {		// only tag=2
		gc_type::npos_vec_type& nv=_gc.bnodes_npos();
		for (size_t k=0; k<nv.size(); k++)
			_set_grad(nv[k].x(),nv[k].y(),true);
	}
}
*/

void TFJSphereAdsorp::prepare_eqn_gama(bool rhs_only)
{
	size_t n=sn();
	_gcoefA.resize(n);	_gcoefB.resize(n);
	_gcoefC.resize(n);	_gcoefR.resize(n);
	params_type& p=_params;
	if (!rhs_only)
	{
		float_t ds=_surf.arcs().seg_len(0);
		float_t t1=p.Dsb*p.hob/ds/ds;
		float_t t2=1.0/p.dt;
		reset_vector(_gcoefA,-t1);
		reset_vector(_gcoefB,t2+2.0*t1);
		reset_vector(_gcoefC,-t1);
		// correct the end nodes
		_gcoefA[0]=0;	_gcoefB[0]=t2+t1;
		int k=n-1;
		_gcoefC[k]=0;	_gcoefB[k]=t2+t1;
	}
	for (size_t i=0; i<n; i++)
		_gcoefR[i]=_pgama[i]/p.dt+_cfluxsub[i];
}

void TFJSphereAdsorp::solve_gama()
{
	_gama_solver.solve();
}

void TFJSphereAdsorp::_dampen_vec(const vector_f& ov, 
							vector_f& nv,float_t alpha)
{
	for (size_t i=0; i<ov.size(); i++)
		nv[i]=alpha*nv[i]+(1-alpha)*ov[i]; 
}

void TFJSphereAdsorp::update_sublayer_conc()
{
	params_type &p=_params;
	float_t k=p.x_eq/(1.0-p.x_eq);	// langmiur model
	float_t tmp=1.0/(1.0-p.x_eq);
	_csub.resize(sn());

	for (size_t i=0; i<sn(); i++)
		_csub[i]=(_cfluxsub[i]/p.lamda+_gama[i])/
				(tmp-k*_gama[i]);
//		_csub[i]=(tmp*_cfluxsub[i]+_gama[i])/
//				(p.k*(1.0/p.x_eq-_gama[i]));	// gama
}

void TFJSphereAdsorp::_set_ghosts()
{
	size_t ii,jj;
	gc_type::npos_vec_type& nv=_gc.ghosts_npos();
	for (size_t i=0; i<nv.size(); i++)
	{
		ii=nv[i].x();	jj=nv[i].y();
		_conc(ii,jj)=_cghost[i];
	}
}

void TFJSphereAdsorp::interp_mirror_csub()
{
	_mcsub.resize(gn());
	vector_f& av=_gc.mirrors_arc();
	_sc_interp.spline();
	for (size_t i=0; i<gn(); i++)
		_mcsub[i]=_sc_interp.splint(av[i]);
}

float_t TFJSphereAdsorp::_dist(float_t x1, float_t y1,
								float_t x2, float_t y2)
{
	float_t dx=x1-x2, dy=y1-y2;
	return sqrt(dx*dx+dy*dy);
}

void TFJSphereAdsorp::interp_cghost()
{
	interp_mirror_csub();
	// prepare
	_gc.set_env(&_genv);
	_gc.fill_env_nodes(_conc,_csub);
	_gc.gen_interp_coef(_g_src,_gcoef);
	// interp
	gc_type::pos_vec_type& mp=_gc.mirrors_pos();
	gc_type::pos_vec_type& gp=_gc.ghosts_pos();
	float dist;
	_mn1v.resize(gn());
	_mcfluxsub.resize(gn());
	vector_f tmp=_mcfluxsub;
	for (size_t i=0; i<gn(); i++)
	{
		_mn1v[i]=_gc.interp_xy(i,_mn1[i].x(),_mn1[i].y(),_gcoef);
		_mcfluxsub[i]=(_mn1v[i]-_mcsub[i])*_params.npl;
	}
//	_dampen_vec(tmp,_mcfluxsub,_params.alpha);
	for (size_t i=0; i<gn(); i++)
	{
		dist=_dist(gp[i].x(),gp[i].y(),mp[i].x(),mp[i].y());
		_cghost[i]=_mcsub[i]-_mcfluxsub[i]*dist;
	}
/*	interp using Ferziger way	
	_mgconc.resize(gn());
	gc_type::pos_vec_type& mp=_gc.mirrors_pos();
	for (size_t i=0; i<gn(); i++)
	{
		_mgconc[i]=_gc.interp_xy(i,_mgpos[i].x(),_mgpos[i].y(),_gcoef);
		_cghost[i]=2.0*_mcsub[i]-_mgconc[i];
	}
*/
	_set_ghosts();
}

/*
void TFJSphereAdsorp::interp_cgrad(const matrix_f& m, vector_f& v)
{
	// prepare
	_gc.set_env(&_senv);
	_gc.fill_env_nodes(m,_csub);	// latter doesn't matter
	_gc.gen_interp_coef(_s_src,_scoef);

	// interp
	size_t ii,x,y;
	v.resize(sn());
	for (size_t i=0; i<sn(); i++)
		v[i]=_gc.interp_xy(i,_surf.x(i),_surf.y(i),_scoef);
}

void TFJSphereAdsorp::_update_surf_bc(vector_f& v)
{
	v[0]=v[1];
	v[sn()-1]=v[sn()-2];
}
*/

void TFJSphereAdsorp::interp_sublayer_cflux()
{
	
	_gc.set_env(&_senv);
	_gc.fill_env_nodes(_conc,_csub);
	_gc.gen_interp_coef(_s_src,_scoef);
//	print_mat("_scoef",_scoef);
	_cn1v.resize(sn());
	_cfluxsub.resize(sn());
	vector_f tmp=_cfluxsub;
	for (size_t i=0; i<sn(); i++)
		_cn1v[i]=_gc.interp_xy(i,_cn1[i].x(),
							_cn1[i].y(),_scoef);
	for (size_t i=0; i<sn(); i++)
		_cfluxsub[i]=(_cn1v[i]-_csub[i])*_params.npl;
	_dampen_vec(tmp,_cfluxsub,_params.alpha);
/*
	if (_scoef.size2()!=6) {
		// linear interpolation scheme
		for (size_t i=0; i<sn(); i++)
			_cfluxsub[i]=-_ysp[i]*_scoef(i,1)+_xsp[i]*_scoef(i,2);
	} else {
		float_t xp,yp;
		for (size_t i=0; i<sn(); i++)
		{
			xp=_scoef(i,1)+_scoef(i,3)*_x[i]*2.0+_scoef(i,4)*_y[i];
			yp=_scoef(i,2)+_scoef(i,4)*_x[i]+_scoef(i,5)*_y[i]*2.0;
			cout << i << " " << xp << " " << yp << endl;
			_cfluxsub[i]=-_ysp[i]*xp+_xsp[i]*yp;
		}
	}
*/
/* using curve node
	interp_cgrad(_cgrad_x,_cfluxsub_x);
	interp_cgrad(_cgrad_y,_cfluxsub_y);
	// assemble cflux, (-n) dot (-DgradC)
	_cfluxsub.resize(sn());
	for (size_t i=0; i<sn(); i++)
		_cfluxsub[i]=-_ysp[i]*_cfluxsub_x[i]+_xsp[i]*_cfluxsub_y[i];
	_update_surf_bc(_cfluxsub);
*/
/* using mirror nodes
	interp_cgrad(_cgrad_x,_mcfluxsub_x);
	interp_cgrad(_cgrad_y,_mcfluxsub_y);
	_mcfluxsub.resize(gn());
	for (size_t i=0;i<gn(); i++)
		_mcfluxsub[i]=-_mysp[i]*_mcfluxsub_x[i]+_mxsp[i]*_mcfluxsub_y[i];
*/
}

void TFJSphereAdsorp::prepare_eqn_conc(bool ghost_only)
{
	size_t m=xn(), n=yn();
	_cap.resize(m,n);		_csu.resize(m,n);
	_caw.resize(m,n);		_cae.resize(m,n);
	_cas.resize(m,n);		_can.resize(m,n);

	params_type &p=_params;
	if (!ghost_only)
	{
	// generate coefficient for bulk
	float_t tmp1=p.hob/dx()/dx(),
			tmp2=p.hob/dy()/dy();;
	reset_matrix(_cap,1.0/p.dt+tmp1*2+tmp2*2);
	reset_matrix(_cas,-tmp2);
	reset_matrix(_can,-tmp2);

	// correction for cylindrical coordinates
	float rp,re,rw;
	for (size_t i=0; i<m; i++)
	{
		rp=_axis.x().seg_center_value(i);
		re=_axis.x().seg_upper_value(i);
		rw=_axis.x().seg_lower_value(i);
		for (size_t j=0; j<n; j++)
		{
			_cae(i,j)=-re/rp*tmp1;
			_caw(i,j)=-rw/rp*tmp1;
			_csu(i,j)=-_pconc(i,j)/p.dt;
		}	
	}
	
	// generate coefficient for regular b.c
	size_t i,j;
	// left and right
	for (j=0; j<n; j++)
	{
		i=0;	_cap(i,j)+=_caw(i,j);	_caw(i,j)=0.0;		
		i=m-1;	_cap(i,j)+=_cae(i,j);	_cae(i,j)=0.0;
	}
	// top and bottom
	for (i=0; i<m; i++)
	{
		j=0;	_cap(i,j)+=_cas(i,j);	_cas(i,j)=0.0;
		j=n-1;	_cap(i,j)+=_can(i,j);;	_can(i,j)=0.0;
	}
	
	}	// end of ghost_only
	
	// generate coefficient for ghost b.c
	gc_type::npos_vec_type& pv=_gc.ghosts_npos();
	size_t ii,jj;
	size_t k=0;
	while (k<pv.size()) {
		ii=pv[k].x(); jj=pv[k].y();
		_caw(ii,jj)=_cae(ii,jj)=_cas(ii,jj)=_can(ii,jj)=0.0;
		_cap(ii,jj)=1.0;
		_csu(ii,jj)=-_conc(ii,jj);	
		k++;
	}
}

void TFJSphereAdsorp::solve_conc()
{
#ifdef _Petsc_KSP
	_cksp.initialize();
	_cksp.solve();
#endif
}

void TFJSphereAdsorp::reset()
{
	matrix_n& m=_gc.get_umat();
	for (size_t i=0; i<xn(); i++)
		for (size_t j=0; j<yn(); j++)
			if (m(i,j)!=0) _conc(i,j)=1.0;
			else _conc(i,j)=0.0;
//	reset_matrix(_conc,1.0);
	reset_vector(_csub,1.0);
	reset_vector(_gama,0.0);
	reset_vector(_cfluxsub,0.0);
	_tc=0.0;
}

float_t TFJSphereAdsorp::_get_cghost_err(const vector_f& vb)
{
	float_t err=0,tmp;
	for (size_t i=0; i<gn(); i++)
	{
		tmp=fabs(_cghost[i]-vb[i]);
		if (tmp>err) err=tmp;
	}
	return err;
}

float_t TFJSphereAdsorp::_get_csub_err(const vector_f& vb)
{
	float_t err=0,tmp;
	for (size_t i=0; i<gn(); i++)
	{
		tmp=fabs(_csub[i]-vb[i]);
		if (tmp>err) err=tmp;
	}
	return err;
}

/*
void TFJSphereAdsorp:: _restrict_vec(vector_f& v)
{
	for (size_t i=0; i< v.size(); i++)
		if (v[i]<0) { 
			cout << "restrict " << v[i] << " to be zero" << endl;
			v[i]=0;
		}
}
*/

void TFJSphereAdsorp::step()
{
	// initialize for the run
	_pgama=_gama;
	_pconc=_conc;
//	_set_ghosts();

	prepare_eqn_conc();
//	update_grad_conc();
	interp_sublayer_cflux();
	prepare_eqn_gama();
/*	
	print_var("_gama",_gama);
	print_var("_cn1v",_cn1v);
	print_var("_csub",_csub);
	print_var("_cfluxsub",_cfluxsub);
*/
//	return;

	_cghost.resize(gn());
	reset_vector(_cghost,0.0);
	vector_f _csub_in,_cghost_out;
	_inner_iters.resize(0);
//	size_t out_iters=0;
	float_t err;
	do
	{
		_cghost_out=_cghost;
		size_t in_iters=0;
		do
		{
			_csub_in=_csub;
			interp_sublayer_cflux();
			prepare_eqn_gama(true);	// prepare rhs only
			solve_gama();
			update_sublayer_conc();
			
			if (_params.debug==1) {
			print_var("in_iters",in_iters);
			print_var("_cn1v",_cn1v);
			print_var("_cfluxsub",_cfluxsub);
			print_var("_gama",_gama);
			print_var("_csub",_csub);
			}
			
			in_iters++;
			err=_get_csub_err(_csub_in);
//		} while (in_iters<3);
		} while	((in_iters==1) || (err>_params.iter_acc));
/*
			print_var("_cfluxsub",_cfluxsub);
			print_var("_gama",_gama);
			print_var("_csub",_csub);
*/
		push_back(_inner_iters,(int)in_iters);
		interp_cghost();
		if (_params.debug>1) {
		print_var("_mn1v",_mn1v);
		print_var("_mcsub",_mcsub);
		print_var("_mcfluxsub",_mcfluxsub);
//		print_var("_mgconc",_mgconc);
		print_var("_cghosts",_cghost);
//		print_mat("_conc",_conc,vmpfText | vmpfColMajor | vmpfColReverse);
		}
//		break;
		prepare_eqn_conc(true);		// prepare ghosts only
		solve_conc();
		if (_params.debug==3) {
			save_mat("_conc.txt",_conc);
			break;
		}
		err=_get_cghost_err(_cghost_out);
		if (_params.debug)
			cout << "o " << _inner_iters.size() << ":" 
				<< in_iters << " " << err << endl;
//		out_iters++;
	} while (err>_params.iter_acc);	// check _tconc,_conc
	_tc+=_params.dt;
}

}	// end of namespace
