#include "fjapp_bubbleNSGD.h"
#include "fjlib_vecmat_print.h"

namespace fjlib {

void TFJBubbleNSGD_Solver::step()
{
	// prepare solve mass
	get_interface_uv();
	_pgama=_gama;
//	print_vec("gama",_gama);
	_pconc=_conc;
	prepare_gc();
	// solve mass
	solve_coupled_mass();
#ifdef _Petsc_KSP
	reset_cksp();
#endif
	// post-solver
	update_surface_tension();
	update_tracked_mass();

	TFJBubbleNS_Solver::step();
}

void TFJBubbleNSGD_Solver::initialize()
{
	TFJBubbleNSG_Solver::initialize();
	_fluxInterp.set_data(_surf.curve().arcs().get_vector(),
						&_cfluxsub);
	// initialize ghost cell method
	initialize_gc();
	// intialize bulk property //	
	_axis.set_Stag(sptO);
	_conc.resize(xc(),yc());
	// initial parallel
#ifdef _Petsc_KSP
	initialize_cksp();
#endif
}

void TFJBubbleNSGD_Solver::reset()
{
	TFJBubbleNSG_Solver::reset();
	if (_tc==0)
		reset_matrix(_conc,1.0);
	reset_vector(_csub,1.0);
	reset_vector(_cfluxsub,0.0);
}

void TFJBubbleNSGD_Solver::initialize_gc()
{
	// setup interpolation gama for mirror
	_sc_interp.set_data(_surf.curve().arcs().get_vector(),&_csub);
//	_sc_interp.set_bc(at::csbtClamped,at::csbtClamped,0,0);
	// setup _cn1 _mn1 src vector 
	set_vector(_cn1src,0,1,2);
	set_vector(_mn1src,0,1,2);
}

void TFJBubbleNSGD_Solver::prepare_gc()
{
	// prepare cn1
	_cn1.resize(sc());
	vector_f& x=_surf.x(), 		y=_surf.y();
	vector_f& xsp=_surf.xp(1),	ysp=_surf.yp(1);
	for (size_t i=0; i<sc()-1; i++)
	{
		_cn1[i].redim(2);
		_cn1[i].x()=x[i]-ysp[i]/_pms.snpl;
		_cn1[i].y()=fabs(y[i]+xsp[i]/_pms.snpl);	// fabs
	}
	// special treatment for the pinned point
	size_t k=sc()-1;
	_cn1[k].redim(2);
	_cn1[k].x()=x[k]+1.0/_pms.snpl;
	_cn1[k].y()=y[k];
/*	
	cout << "_x,_y,_cn1x,_cn1y:" << endl;
	for (size_t i=0; i<sc(); i++)
		cout << i 
			<< " " << _surf.x()[i] << "," << _surf.y()[i] 
			<< " " << _cn1[i].x() << "," << _cn1[i].y() << endl;
*/	
	curve_type& cv=_surf.curve();
	// setup gc
	_vof.fill_unknowns(_pre_umat);
	_gc.set_data(&_axis.axis(sptO),&cv,_pre_umat);
	_gc.gc_gen();
//	save_mat("umat.txt",_gc.get_umat());
	
	//  prepare mn1
	_mxsp.resize(gc()); _mysp.resize(gc());
	_mn1.resize(gc());
	if (gc()>0) {
	vector_f tmp;
	for (size_t i=0; i<gc(); i++)
	{
		gc_type::pos_type a=_gc.mirrors_pos()[i];
		gc_type::pos_type b=_gc.ghosts_pos()[i];
		
		float_t ss=_gc.mirrors_arc()[i];
		if (ss>=0) {
		cv.interps().x().splint(ss,&tmp);	_mxsp[i]=tmp[1];
		cv.interps().y().splint(ss,&tmp);	_mysp[i]=tmp[1];
		
		_mn1[i].redim(2);
		_mn1[i].x()=a.x()-_mysp[i]/_pms.snpl;
		_mn1[i].y()=a.y()+_mxsp[i]/_pms.snpl;
		}
	}
	}	// end of if(gc()>0)
/*
	cout << "_mn1x,_mn1y:" << endl;
	for (size_t i=0; i<gc(); i++)
		cout << i 
			<< " " << _mn1[i].x() << "," << _mn1[i].y() << endl;
*/
	// prepare _cn1, _mn1 env
	_gc.set_env(&_cn1env);
	_gc.gen_env_info(_cn1,3,0);		// for flux
	if (gc()>0) {
	_gc.set_env(&_mn1env);
	_gc.gen_env_info(_mn1,3,0);		// for cg
	// make conc look nice, remove void
	// optimize later
	matrix_n& m=_gc.get_umat();
	for (size_t i=0; i<m.size1(); i++)
		for (size_t j=0; j<m.size2(); j++)
			if (m(i,j)==0) _conc(i,j)=0.0;
	}
}

void TFJBubbleNSGD_Solver::prepare_eqn_gama_rhs()
{ 
//	if (!_diffusion_switch)
	if (_dpms.rate_limit_step<2)
		TFJBubbleNSG_Solver::prepare_eqn_gama_rhs();
	else {
		for (size_t i=0; i<(int)sc()-1; i++)
			_coefR[i]+=_cfluxsub[i]/_dpms.Peb/_dpms.hob;
	}
}

void TFJBubbleNSGD_Solver::_dampen_vec(const vector_f& ov, 
							vector_f& nv,float_t alpha)
{
	for (size_t i=0; i<ov.size(); i++)
		nv[i]=alpha*nv[i]+(1-alpha)*ov[i]; 
}

void TFJBubbleNSGD_Solver::guess_sublayer_conc()
{
	_csub.resize(sc());
	float_t k=c0oa();	
	if (lamda()>1)	// diffusion controlled
		for (size_t i=0; i<sc(); i++)
			_csub[i]=_gama[i]/k/(1.0/_gpms.x0-_gama[i]);
	else			// adsorption controlled or mixed
		reset_vector(_csub,1.0);
}

void TFJBubbleNSGD_Solver::update_sublayer_conc()
{
	float_t k=c0oa();	
	float_t lam=_gpms.Bi*_dpms.Peb*_dpms.hob;
	_csub.resize(sc());
	switch (_dpms.rate_limit_step) {
		case rlsDiffusion:
			for (size_t i=0; i<sc(); i++)
				_csub[i]=_gama[i]/k/(1.0/_gpms.x0-_gama[i]);
			break;
		case rlsMixed:
			for (size_t i=0; i<sc(); i++)
				_csub[i]=(_cfluxsub[i]/lam+_gama[i])/k/
					(1.0/_gpms.x0-_gama[i]);
			break;
		default:
			throw "please select correct rate limit step type";
	}
}

void TFJBubbleNSGD_Solver::_set_ghosts()
{
	if (gc()==0) return;
	size_t ii,jj;
	gc_type::npos_vec_type& nv=_gc.ghosts_npos();
	for (size_t i=0; i<nv.size(); i++)
	{
		ii=nv[i].x();	jj=nv[i].y();
		_conc(ii,jj)=_cghost[i];
	}
}

void TFJBubbleNSGD_Solver::interp_mirror_csub()
{
	_mcsub.resize(gc());
	vector_f& av=_gc.mirrors_arc();
	_sc_interp.spline();
	for (size_t i=0; i<gc(); i++)
		if (av[i]>=0)
			_mcsub[i]=_sc_interp.splint(av[i]);
}

float_t TFJBubbleNSGD_Solver::_dist(float_t x1, float_t y1,
								float_t x2, float_t y2)
{
	float_t dx=x1-x2, dy=y1-y2;
	return sqrt(dx*dx+dy*dy);
}

void TFJBubbleNSGD_Solver::interp_cghost()
{
	if (gc()==0) return;
	interp_mirror_csub();
	// prepare
	_gc.set_env(&_mn1env);
	_gc.fill_env_nodes(_conc,_csub);
	_gc.gen_interp_coef(_mn1src,_mn1coef);
	// interp
	gc_type::pos_vec_type& mp=_gc.mirrors_pos();
	gc_type::pos_vec_type& gp=_gc.ghosts_pos();
	float dist;
	_mn1v.resize(gc());
	_mcfluxsub.resize(gc());
	vector_f tmp=_mcfluxsub;
	for (size_t i=0; i<gc(); i++)
	{
		float_t ss=_gc.mirrors_arc()[i];
		if (ss>=0) {
			_mn1v[i]=_gc.interp_xy(i,_mn1[i].x(),
									_mn1[i].y(),_mn1coef);
			_mcfluxsub[i]=(_mn1v[i]-_mcsub[i])*_pms.snpl;
			dist=_dist(gp[i].x(),gp[i].y(),mp[i].x(),mp[i].y());
			_cghost[i]=_mcsub[i]-_mcfluxsub[i]*dist;
		}
		else if (i!=0)
			throw "missing mirror node in fjapp_bubbleNSGD.cpp";
	}
	size_t gi=0;	// locate bottom ghost node id
	if (_gc.mirrors_arc()[gi]<0)
	{
		cout << "fixed ghost node interpolation at bottom." << endl;
		// take the pinned surf node value and flux
		_mcfluxsub[gi]=_mcfluxsub[gi+1];
		size_t k=sc()-1;
		_mn1v[gi]=_csub[k];	// use this variable to store
		dist=_dist(gp[gi].x(),gp[gi].y(),
					_surf.x()[k],_surf.y()[k]);
		_cghost[gi]=_mn1v[gi]-_mcfluxsub[gi]*dist;
		// need to correct the intersect angle later for accuracy
	}
	_set_ghosts();
}

void TFJBubbleNSGD_Solver::interp_sublayer_cflux()
{
	_gc.set_env(&_cn1env);
	_gc.fill_env_nodes(_conc,_csub);
	_gc.gen_interp_coef(_cn1src,_cn1coef);
	_cn1v.resize(sc());
	_cfluxsub.resize(sc());
	vector_f tmp=_cfluxsub;
	for (size_t i=0; i<sc(); i++)
		_cn1v[i]=_gc.interp_xy(i,_cn1[i].x(),
							_cn1[i].y(),_cn1coef);
	for (size_t i=0; i<sc(); i++)
		_cfluxsub[i]=(_cn1v[i]-_csub[i])*_pms.snpl;
	// special treatment for the pinned node
	size_t k=sc()-1;
	_cfluxsub[k]*=-_surf.yp(1)[k];
	_dampen_vec(tmp,_cfluxsub,_dpms.alpha);
}

void TFJBubbleNSGD_Solver::prepare_eqn_conc(bool ghost_only)
{
	_axis.set_Stag(sptO);
	size_t m=xc(), n=yc();
	_cap.resize(m,n);		_csu.resize(m,n);
	_caw.resize(m,n);		_cae.resize(m,n);
	_cas.resize(m,n);		_can.resize(m,n);

	params_diff_type &p=_dpms;
	if (!ghost_only)
	{
	float_t dr=dx(0), dz=dy(0);
	// generate coefficient for bulk
	float_t tmp1=1.0/p.Peb/dr/dr,
			tmp2=1.0/p.Peb/dz/dz;
	reset_matrix(_cap,1.0/_pms.dt+tmp1*2+tmp2*2);
	reset_matrix(_cas,-tmp2);
	reset_matrix(_can,-tmp2);

	// correction for cylindrical coordinates
	float rp,re,rw;
	float ue,uw,vn,vs;
	axis_type ax=_axis.axis(sptO);
	for (size_t i=0; i<m; i++)
	{
		rp=ax.x().seg_center_value(i);
		re=ax.x().seg_upper_value(i);
		rw=ax.x().seg_lower_value(i);
		for (size_t j=0; j<n; j++)
		{
			ue=_u(i+1,j+1);	uw=_u(i,j+1);
			vn=_v(i+1,j+1); vs=_v(i+1,j);
			_cap(i,j)+=(re*ue-rw*uw)/dr/rp/2.0+
						(vn-vs)/dz/2.0;
			_cae(i,j)=-re/rp*tmp1+re/rp/dr/2.0*ue;
			_caw(i,j)=-rw/rp*tmp1-rw/rp/dr/2.0*uw;
			_can(i,j)+=1.0/dz/2.0*vn;
			_cas(i,j)+=-1.0/dz/2.0*vs;
			_csu(i,j)=-_pconc(i,j)/_pms.dt;
		/*
			if (_gc.get_umat()(i,j)!=0)  // check unknown
				if ((_cap(i,j)<0) || (_caw(i,j)>0) ||
					(_cae(i,j)>0) || (_can(i,j)>0) ||
					(_cas(i,j)>0))
				{
					cout << "_cap,_caw,_cae,_cas,_can at " 
						<< i << "," << j << endl;
					cout <<_cap(i,j) << " " << _caw(i,j) << " "
						<< _cae(i,j) << " " << _cas(i,j) << " "
						<< _can(i,j);
					throw "coef for conc are not bounded";
				}
		*/
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

float_t TFJBubbleNSGD_Solver::_get_cghost_err(const vector_f& vb)
{
	float_t err=0,tmp;
	for (size_t i=0; i<gc(); i++)
	{
		tmp=fabs(_cghost[i]-vb[i]);
		if (tmp>err) err=tmp;
	}
	return err;
}

float_t TFJBubbleNSGD_Solver::_get_csub_err(const vector_f& vb)
{
	float_t err=0,tmp;
	for (size_t i=0; i<gc(); i++)
	{
		tmp=fabs(_csub[i]-vb[i]);
		if (tmp>err) err=tmp;
	}
	return err;
}

void TFJBubbleNSGD_Solver::solve_coupled_mass()
{
	_csub.resize(sc());
	_cfluxsub.resize(sc());
	prepare_eqn_conc();
//	interp_sublayer_cflux();
//	prepare_eqn_gama();
	
	_cghost.resize(gc());
	reset_vector(_cghost,1.0);
	vector_f _csub_in,_cghost_out;
	_inner_iters.resize(0);
	float_t err;
//	reset_vector(_cfluxsub,0.0);
//	guess_sublayer_conc()
	reset_vector(_csub,1.0);
	do
	{
		_cghost_out=_cghost;
		size_t in_iters=0;
		do
		{
			_csub_in=_csub;
			interp_sublayer_cflux();
//			prepare_eqn_gama(true);	// prepare rhs only
			solve_surface_mass();
			update_sublayer_conc();
			in_iters++;
			err=_get_csub_err(_csub_in);
		} while	((in_iters==1) || (err>_dpms.iter_acc));
		push_back(_inner_iters,(int)in_iters);
		/*
			print_vec("iters",_inner_iters);
			print_vec("cn1v",_cn1v);
			print_vec("cflux",_cfluxsub);
			print_vec("gama",_gama);
			print_vec("csub",_csub);
		*/
		interp_cghost();
//			print_vec("cghost",_cghost);
		prepare_eqn_conc(true);		// prepare ghosts only
		solve_conc();
		err=_get_cghost_err(_cghost_out);
	} while (err>_dpms.iter_acc);	// check _tconc,_conc
/*
	cout << '\t' << '\t' << "iters:";
	for (size_t i=0; i<_inner_iters.size(); i++)
		cout << " " << _inner_iters[i];
	cout << endl;
*/
//			save_mat("conc.txt",_conc);
}

float_t TFJBubbleNSGD_Solver::get_bulk_mass()
{
	_axis.set_Stag(sptO);
	size_t m=xc(), n=yc();
	matrix_n& um=_gc.get_umat();
	float_t _mass=0.0,xp;
	for (size_t i=0; i<m; i++)
		for (size_t j=0; j<n; j++)
			if (um(i,j)>0) {
				_mass+=(sqr(xe(i))-sqr(xw(i)))*dy(j)*_conc(i,j);
			}
	return _mass*M_PI;
}

float_t TFJBubbleNSGD_Solver::get_bulk_mass_err()
{
	float_t mt=get_bulk_mass()+get_curve_mass();
	float_t m0=M_PI*sqr(_pms.xmax)*_pms.ymax;
	if (_pms.flat) m0+=M_PI;
	else m0+=2.0*M_PI;
	return (mt-m0)/m0;
}

float_t TFJBubbleNSGD_Solver::get_absorbed_mass()
{
	float_t _mass=0.0,ds;
	TFJPolynomial p,g,pl;
	g.set_degree(1);
	for (size_t i=0; i<_surf.curve().seg_count(); i++)
	{
		ds=_surf.curve().arcs().seg_len(i);
		g.set_coef(_cfluxsub[i+1]-_cfluxsub[i],_cfluxsub[i]);
		p=_surf.curve().interps().x().get_seg_poly(i);

        pl=(p*g)<<1;
		_mass+=pl.eval_diff(0.0,1.0)*ds;
	}
	return 2.0*M_PI*_mass/_dpms.Peb/_dpms.hob*_pms.dt;
}

float_t TFJBubbleNSGD_Solver::lamda()
{
	float_t ld;
	switch (_dpms.rate_limit_step) {
		case rlsDiffusion:
			ld=1000;
			break;
		case rlsAdsorption:
			ld=0.001;
			break;
		case rlsMixed:
			ld=_gpms.Bi/(1.0-_gpms.x0)*_dpms.Peb*_dpms.hob;
			break;
		default:
			ld=1;
			break;
	}
	return ld;
}

#ifdef _Petsc_KSP
void TFJBubbleNSGD_Solver::initialize_cksp()
{
	_cksp.set_data(&_cap,&_caw,&_cae,
					&_cas,&_can,&_csu,&_conc);
	_cksp.set_mask(&_gc.get_umat());
	_cksp.set_tolerance(_dpms.conc_acc);
}

void TFJBubbleNSGD_Solver::reset_cksp()
{
	_cksp.reset();
}

void TFJBubbleNSGD_Solver::solve_conc()
{
	_cksp.initialize();
	_cksp.solve();
}
#endif

}	// end of namespace
