#include "fjapp_bubbleNSG.h"

namespace fjlib {

void TFJBubbleNSG_Solver::step()
{
	// prepare solve mass
	get_interface_uv();
	_pgama=_gama;
	// solve mass
	solve_surface_mass();
	// post-solver
	update_surface_tension();
	update_tracked_mass();

	TFJBubbleNS_Solver::step();
}

void TFJBubbleNSG_Solver::initialize_interface()
{
	TFJBubbleNS_Solver::initialize_interface();
	// set interpolation for gama and sigma
	_gInterp.set_data(_surf.curve().arcs().get_vector(),&_gama);
	_sigInterp.set_data(_surf.curve().arcs().get_vector(),&_sigma);
	// set surfactant data
	_surfEq.set_params(_gpms.Surfc);
	_surfEq.set_Gama(_gpms.x0);
	_surfEq.calc_Sigma();
	_surfEq.calc_C();
	_surfDyn.set_Eq(&_surfEq);
	_max_gama=_surfEq.solve_MaxGama()/_surfEq.get_Gama();
	_gama_solver.set_data(&_coefB,&_coefA,&_coefC,&_coefR,&_gama);

	_vof.set_STInterp(&_sigInterp);
}

void TFJBubbleNSG_Solver::reset_interface()
{
	TFJBubbleNS_Solver::reset_interface();
	if (_tc==0)
	{
		size_t cn=_surf.curve().pt_count();
		_gama.resize(cn);
		reset_vector(_gama,1.0);
	}
//	else should be loaded manually

	// initial tracked mass
	_tracked_mass=get_curve_mass();
}

void TFJBubbleNSG_Solver::remesh_interface()
{
//	_surf.remesh(_axis.axis(sptO));
	_surf.remesh_uniform();
	remesh_gama();
	update_surface_tension();
	_surf.finalize_mesh();
	_surf.spline();
}

void TFJBubbleNSG_Solver::update_surface_tension()
{
	_sigma.resize(_gama.size());
	for (size_t i=0; i<_gama.size(); i++)
	{
		// make sure _gama is bounded below maxGama
		if (_gama[i]>_max_gama) _gama[i]=_max_gama;
		_surfDyn.set_Gama(_gama[i]);
		_surfDyn.calc_Sigma();
		_sigma[i]=_surfDyn.get_Sigma()*_surfEq.get_Sigma();	// scaled with sigma0
	}
}

void  TFJBubbleNSG_Solver::remesh_gama()
{
	vector_f& ns=_surf.new_mesh_arc();
	vector_f ng(ns.size());
	for (size_t i=0; i<ns.size(); i++)
		ng[i]=_gInterp.splint(ns[i]);
	_gama=ng;
}

/*
void TFJBubbleNSG_Solver::solve_interface()
{
	get_interface_uv();
	move_interface();
	prepare_interface();
	get_interface_uv();
	solve_surface_mass();
}
*/

void TFJBubbleNSG_Solver::prepare_eqn_gama2()
{
	_surf.calc_curvature(_cvt);

	size_t cn=_surf.curve().pt_count();
    _coefA.resize(cn); _coefB.resize(cn);
	_coefC.resize(cn); _coefR.resize(cn);

	// higher order scheme preparation
	vector_f gupper,gdown;
	gupper.resize(cn); gdown.resize(cn);
	float_t ghat;
	for (size_t i=0; i<cn-2; i++)
	{
		ghat=gama_HatD2N(_gama[i+1],_gama[i],_gama[i+2]);
		gupper[i+1]=gama_HatN2D(gama_H2F(ghat),_gama[i],_gama[i+2]);
		ghat=gama_HatD2N(_gama[i+1],_gama[i+2],_gama[i]);
		gdown[i]=gama_HatN2D(gama_H2F(ghat),_gama[i+2],_gama[i]);
	}

	float_t dsp,dse,dsw,sp,spe,spw,gh,gf;
	float_t xe,xw,xp,use,usw,de,dw,fe,fw;
	float_t ppe=_gpms.Pe;
	// domain
	for (size_t i=1; i<cn-1; i++)
	{	
		// calculating seg length
		dse=_surf.curve().arcs().seg_len(i);
		dsw=_surf.curve().arcs().seg_len(i-1);
		dsp=(dse+dsw)/2.0;
		sp=_surf.curve().arcs().value_at(i);
		spe=sp+dse/2.0;
		spw=sp-dsw/2.0;

		xp=_surf.curve().x(i);
		xe=_surf.curve().interp(0).splint(spe);
		xw=_surf.curve().interp(0).splint(spw);

		use=(_ut[i]+_ut[i+1])/2.0;
		usw=(_ut[i]+_ut[i-1])/2.0;

		fe=xe*use/dsp/xp;		fw=xw*usw/dsp/xp;
		de=xe/dse/dsp/ppe/xp;	dw=xw/dsw/dsp/ppe/xp;
//		de=0.0; dw=0.0;  // no diffusion test

		_coefA[i]=-dw-std::max(fw,0.0);
		_coefC[i]=-de-std::max(-fe,0.0);
		_coefR[i]=_pgama[i]/_pms.dt*2.0;
		_coefB[i]=-(_coefA[i]+_coefC[i])+(fe-fw);
		_coefB[i]+=_un[i]*_cvt[i]+2.0/_pms.dt;

		// include higher order correction
		_coefR[i]-=std::max(fe,0.0)*(gupper[i]-_gama[i]);
		if (i<cn-2)
			_coefR[i]+=std::max(-fe,0.0)*(gdown[i]-_gama[i+1]);
		if (i>1)
			_coefR[i]+=std::max(fw,0.0)*(gupper[i-1]-_gama[i-1]);
		_coefR[i]-=std::max(-fw,0.0)*(gdown[i-1]-_gama[i]);

		// crank-nicolson scheme
		_coefR[i]-=(-dw-fw/2.0)*_pgama[i-1]+
					(-de+fe/2.0)*_pgama[i+1]+
					(de+dw+(fe-fw)/2.0+_un[i]*_cvt[i])*_pgama[i];
	}
	//s0
	_coefA[0]=0.0;
	_coefB[0]=2.0/_pms.dt+_un[0]*_cvt[0];
	_coefC[0]=0.0;
	_coefR[0]=_pgama[0]/_pms.dt*2.0-_un[0]*_cvt[0]*_pgama[0];
	//send
	_coefA[cn-1]=-1.0;
	_coefB[cn-1]=1.0;
	_coefC[cn-1]=0.0;
	_coefR[cn-1]=0.0;

	if (_gpms.Bi>0)
	{
		float_t coa=c0oa(),
				cov=coverage();
		float_t rr=_gpms.Bi*coa/cov;
		float_t bb=_gpms.Bi*(coa+1);
		for (size_t i=0; i<cn-1; i++)
		{
			_coefB[i]+=bb;
			_coefR[i]+=rr*2.0-bb*_pgama[i];
		}
	}
}

void TFJBubbleNSG_Solver::prepare_eqn_gama(bool rhs_only)
{
	if (!rhs_only)
	{
//	vector_f _cvt;
	_surf.calc_curvature(_cvt);

//	size_t cn=_surf.curve().pt_count();
	size_t cn=sc();
    _coefA.resize(cn); _coefB.resize(cn);
	_coefC.resize(cn); _coefR.resize(cn);

	// higher order scheme preparation
	vector_f gupper,gdown;
	gupper.resize(cn); gdown.resize(cn);
	float_t ghat;
	for (size_t i=0; i<cn-2; i++)
	{
		ghat=gama_HatD2N(_gama[i+1],_gama[i],_gama[i+2]);
		gupper[i+1]=gama_HatN2D(gama_H2F(ghat),_gama[i],_gama[i+2]);
		ghat=gama_HatD2N(_gama[i+1],_gama[i+2],_gama[i]);
		gdown[i]=gama_HatN2D(gama_H2F(ghat),_gama[i+2],_gama[i]);
	}

	float_t dsp,dse,dsw,sp,spe,spw,gh,gf;
	float_t xe,xw,xp,use,usw,de,dw,fe,fw;
	float_t ppe=_gpms.Pe;
	// domain
	for (size_t i=1; i<cn-1; i++)
	{	
		// calculating seg length
		dse=_surf.curve().arcs().seg_len(i);
		dsw=_surf.curve().arcs().seg_len(i-1);
		dsp=(dse+dsw)/2.0;
		sp=_surf.curve().arcs().value_at(i);
		spe=sp+dse/2.0;
		spw=sp-dsw/2.0;

		xp=_surf.curve().x(i);
		xe=_surf.curve().interp(0).splint(spe);
		xw=_surf.curve().interp(0).splint(spw);

		use=(_ut[i]+_ut[i+1])/2.0;
		usw=(_ut[i]+_ut[i-1])/2.0;

		fe=xe*use/dsp/xp;		fw=xw*usw/dsp/xp;
		de=xe/dse/dsp/ppe/xp;	dw=xw/dsw/dsp/ppe/xp;
//		de=0.0; dw=0.0;  // no diffusion test

		_coefA[i]=-dw-std::max(fw,0.0);
		_coefC[i]=-de-std::max(-fe,0.0);
		_coefR[i]=_pgama[i]/_pms.dt;
		_coefB[i]=-(_coefA[i]+_coefC[i])+(fe-fw);
		_coefB[i]+=_un[i]*_cvt[i]+1.0/_pms.dt;
		// include higher order correction
		_coefR[i]-=std::max(fe,0.0)*(gupper[i]-_gama[i]);
		if (i<cn-2)
			_coefR[i]+=std::max(-fe,0.0)*(gdown[i]-_gama[i+1]);
		if (i>1)
			_coefR[i]+=std::max(fw,0.0)*(gupper[i-1]-_gama[i-1]);
		_coefR[i]-=std::max(-fw,0.0)*(gdown[i-1]-_gama[i]);
	}
	//s0
	_coefA[0]=0.0;
	_coefB[0]=1.0/_pms.dt+_un[0]*_cvt[0];
	_coefC[0]=0.0;
	_coefR[0]=_pgama[0]/_pms.dt;
	//send
	_coefA[cn-1]=-1.0;
	_coefB[cn-1]=1.0;
	_coefC[cn-1]=0.0;
	_coefR[cn-1]=0.0;
/*
	//send
	use=(_ut[cn-1]+_ut[cn-2])/2.0;
	dse=_surf.curve().arcs().seg_len(cn-2);
	dsp=dse;
	_coefA[cn-1]=0.0;
	_coefB[cn-1]=1.0/_pms.dt-use/dse;
	_coefC[cn-1]=0.0;
	_coefR[cn-1]=_pgama[cn-1]/_pms.dt;
*/
	}
	prepare_eqn_gama_rhs();
}

void TFJBubbleNSG_Solver::prepare_eqn_gama_rhs()
{
	if (_gpms.Bi>0)
	{
		float_t coa=c0oa(),
				cov=coverage();
		for (size_t i=0; i<(int)sc()-1; i++)
		{
			_coefR[i]+=_gpms.Bi*coa/cov;
			_coefB[i]+=_gpms.Bi*(coa+1);
		}
	}
}

float_t	TFJBubbleNSG_Solver::gama_HatD2N(float_t val,
										float_t upstream,
										float_t downstream)
{
	if (upstream==downstream) return upstream;
	else return (val-upstream)/(downstream-upstream);
}

float_t	TFJBubbleNSG_Solver::gama_HatN2D(float_t val,
										float_t upstream,
										float_t downstream)
{
	return val*(downstream-upstream)+upstream;
}

float_t	TFJBubbleNSG_Solver::gama_H2F(float_t v)
{
// unkown somwhere
//
	if ((v>0) && (v<0.2)) return 3.0*v;
	if (v>=0.2) return (v+1)/2.0;
	return v/2;
//
/*
	if ((v>=0) && (v<=1.0/6.0)) return 3.5*v;
	if ((v>1.0/6.0) && (v<=0.5)) return 0.5*(v+1.0);
	if ((v>0.5) && (v<=0.75)) return v+0.25;
	if ((v>=0.75) && (v<=1.0)) return 1.0;
	return v/2.0;
*/
}

void TFJBubbleNSG_Solver::solve_gama()
{
	size_t cn=_surf.curve().pt_count();
/*
	// convert to tridiagonal matrix
	_coefB[0]=_coefA[1]*_coefA[0]-_coefB[0]*_coefC[1];
	_coefC[0]=_coefB[1]*_coefA[0]-_coefC[0]*_coefC[1];
	_coefR[0]=_coefR[1]*_coefA[0]-_coefR[0]*_coefC[1];
	_coefA[0]=0.0;
	_coefA[cn-1]=_coefB[cn-2]*_coefC[cn-1]-_coefA[cn-1]*_coefA[cn-2];
	_coefB[cn-1]=_coefC[cn-2]*_coefC[cn-1]-_coefB[cn-1]*_coefA[cn-2];
	_coefR[cn-1]=_coefR[cn-2]*_coefC[cn-1]-_coefR[cn-1]*_coefA[cn-2];
	_coefC[cn-1]=0.0;
*/
	_gama_solver.solve();
}

void TFJBubbleNSG_Solver::solve_surface_mass()
{
//	_pgama=_gama;
	int iter=0;
	float_t tmp,err;
	do {
		vector_f og=_gama;
		prepare_eqn_gama();
		solve_gama();
		err=0.0;
		for (size_t i=0; i<_gama.size(); i++)
		{
			if (_gama[i]!=0)
				tmp=std::fabs((og[i]-_gama[i])/_gama[i]);
			if (tmp>err) err=tmp;
		}		
		iter++;
	} while (err>_gpms.gacc);
//	update_surface_tension();
//	update_tracked_mass();
}

float_t TFJBubbleNSG_Solver::get_curve_mass()
{
	float_t _mass=0.0,ds,ga,gb;
	TFJPolynomial p,g,pl;
	for (size_t i=0; i<_surf.curve().seg_count(); i++)
	{
		ds=_surf.curve().arcs().seg_len(i);
		g=_gInterp.get_seg_poly(i);
		p=_surf.curve().interps().x().get_seg_poly(i);

        pl=(p*g)<<1;
		_mass+=pl.eval_diff(0.0,1.0)*ds;
	}
	return _mass*2.0*M_PI;
}

float_t TFJBubbleNSG_Solver::get_absorbed_mass()
{
	if (_gpms.Bi<=0) return 0;
	float_t _mass=0.0,ds,ga,gb;
	TFJPolynomial p,g,pl;
	for (size_t i=0; i<_surf.curve().seg_count(); i++)
	{
		ds=_surf.curve().arcs().seg_len(i);
		g=_gInterp.get_seg_poly(i);
		p=_surf.curve().interps().x().get_seg_poly(i);

        pl=(p-p*g)<<1;
		_mass+=pl.eval_diff(0.0,1.0)*ds;
	}
	return 2.0*M_PI*_mass*
				_gpms.Bi/(1.0-coverage())*_pms.dt;
}

}		// end of namespace
