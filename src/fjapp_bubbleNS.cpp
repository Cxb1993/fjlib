#include "fjapp_bubbleNS.h"
#include "fjlib_vecmat_print.h"

namespace fjlib {

void TFJBubbleNS_Interface::apply_spline_bc()
{
	typedef TFJCubicSpline a;
	if (type_id==0) 
	{	
		// set the begin bc to be first derivative fixed
		// set the end bc to be natural, second derivative=0
		_curve.interps().x().set_bc(a::csbtClamped,
								a::csbtNatural,1,0);
		_curve.interps().y().set_bc(a::csbtClamped,
								a::csbtNatural,0,0);
	} else 
	{
		_curve.interps().x().set_bc(a::csbtClamped,
								a::csbtClamped,1,-1);
		_curve.interps().y().set_bc(a::csbtClamped,
								a::csbtClamped,0,0);
	}
}

void TFJBubbleNS_Interface::initialize()
{
	_curve.set_data(&_x,&_y);
	set_type(0);
}

void TFJBubbleNS_Interface::reset()
{
	if (!flat)
	{
		size_t nodes=M_PI_2*npl+1;
		if (nodes<2) nodes=2;
		_x.resize(nodes);	_y.resize(nodes);

		float_t dsita=M_PI_2/(nodes-1);
		for (size_t i=0; i<nodes; i++)
		{
			float_t ac=dsita*i;
			_x[i]=std::sin(ac);
			_y[i]=std::cos(ac);
		}
	}
	else
	{
		size_t nodes=npl;
		if (nodes<2) nodes=2;
		_x.resize(nodes);	_y.resize(nodes);

		float_t dsita=1.0/(nodes-1);
		for (size_t i=0; i<nodes; i++)
		{
			float_t ac=dsita*i;
			_x[i]=ac;
			_y[i]=0;
		}
	}
	spline();
}

void TFJBubbleNS_Interface::spline()
{
	apply_bc();
	_curve.calc_arcs();
	_curve.spline();
	calc_derivs();
}

void TFJBubbleNS_Interface::apply_bc()
{
	size_t cn=_curve.pt_count();
	if (cn<1) return;

	// Fix coords since spline will change it
	_x[0]=0;			// changed recently for diffusion
							// don't know why it's not working 
							// for non-flat case
	if (type_id==0)
		_y[cn-1]=0;
	else
		_x[cn-1]=-1e-5;		// hack for identify multiple drops

	_xp1.resize(cn);	_xp2.resize(cn);
	_yp1.resize(cn);	_yp2.resize(cn);
	_arc_crt.resize(cn); 

	// Correction b.c
	_xp1[0]=1;		_yp1[0]=0;
	if (type_id!=0)
	{ 
		_xp1[cn-1]=-1;	_yp1[cn-1]=0;
	}
}

void TFJBubbleNS_Interface::arc_deriv(const vector_f& v, size_t ord, 
									vector_f& ov)
{
	size_t cn=_curve.pt_count();
	if (v.size()!=cn) throw "\nerror 1 in fjapp_bubbleNS.cpp\n";
	ov.resize(cn);
	switch (ord)	{
		case 0:
			ov=v; break;
		case 1:
			if (cn<2) throw "\nerror 2 in fjapp_bubbleNS.cpp\n";
			for (size_t i=1; i<cn-1; i++)
				ov[i]=(v[i+1]-v[i-1])/
						(_curve.arcs().seg_len(i-1)+
						_curve.arcs().seg_len(i));
			ov[0]=ov[1];
			ov[cn-1]=ov[cn-2];
			break;
		case 2:
			if (cn<3) throw "\nerror 3 in fjapp_bubbleNS.cpp\n";
			float_t de,dw;
			for (size_t i=1; i<cn-1; i++)
			{
				de=_curve.arcs().seg_len(i);
				dw=_curve.arcs().seg_len(i-1);
				ov[i]=((v[i+1]-v[i])/de-
						(v[i]-v[i-1])/dw)/(de+dw)*2.0;
			}
			ov[0]=ov[1];
			ov[cn-1]=ov[cn-2];
			break;
	}
}

void TFJBubbleNS_Interface::calc_derivs()
{
	size_t cn=_curve.pt_count();
	_xp1.resize(cn);	_xp2.resize(cn);
	_yp1.resize(cn);	_yp2.resize(cn);
	_arc_crt.resize(cn);

	vector_f nx,ny;
	TFJSmoothing sm_obj;
	sm_obj.set(2,2,2,true);
	sm_obj.smoothing(_x,nx);
	sm_obj.smoothing(_y,ny);
	arc_deriv(nx,1,_xp1);		
	arc_deriv(ny,1,_yp1);

//	SmoothVector(nx);		SmoothVector(ny);
	arc_deriv(nx,2,_xp2);
	arc_deriv(ny,2,_yp2);

	apply_bc();

	for (size_t i=0; i<cn; i++)
		_arc_crt[i]=std::sqrt(sqr(_xp1[i])+sqr(_yp1[i]));
}

bool TFJBubbleNS_Interface::remesh()
{
	// locate one segment need to be adjusted
	float_t ds=1.0/npl;
	int loc=-1;
	bool tolong;
	for (size_t i=0; i<_curve.seg_count(); i++)
	{
		float_t len=_curve.arcs().seg_len(i);
		if (len>2.4*ds)
		{ loc=i; tolong=true; break;}
//		if (len<0.6*ds)
//		{ loc=i; tolong=false; break;}
	}
	if (loc<0) return false;

	if (!tolong)	// too short, then adjust to the middle
	{
		size_t j=loc+1;
		if (j==_curve.pt_count()-1) j--;
		float_t ss=(_curve.arcs().value_at(j-1)+
					_curve.arcs().value_at(j+1))/2.0;
		curve_type::point_type a=_curve.splint(ss);
		_curve.x(j)=a.x();
		_curve.y(j)=a.y();
	}
	else			// too long, insert one in the middle
	{
		size_t j=loc;
		float_t ss=_curve.arcs().seg_center_value(j);
		curve_type::point_type a=_curve.splint(ss);
		insert(_x,j+1,a.x());
		insert(_y,j+1,a.y());
	}  
	return true;
}

void TFJBubbleNS_Interface::remesh_uniform()
{
	size_t nsegs=_curve.arcs().len()*npl;
	float_t ds=_curve.arcs().len()/nsegs;

	size_t npts=nsegs+1;
	_ns.resize(npts);
	for (size_t i=0; i<npts; i++)
		_ns[i]=ds*i;

/*
	size_t npts=nsegs+1;
	vector_f nx(npts),ny(npts);
	curve_type::point_type a;
	for (size_t i=0; i<npts; i++)
	{
		a=_curve.splint(ds*i);
		nx[i]=a.x(); ny[i]=a.y();
	}
	_x=nx;	_y=ny;
*/
}

void TFJBubbleNS_Interface::remesh(const axis_type& axis)
{
	// make a vector recording mesh density of each node
	size_t n=_curve.pt_count();
	vector_f md(n);
	for (size_t i=0; i<n; i++)
	{
		// find the node location ii,jj
		size_t ii=axis.x().locate_seg(_curve.x(i));
		float_t sx=axis.x().seg_len(ii);
		size_t jj=axis.y().locate_seg(_curve.y(i));
		float_t sy=axis.y().seg_len(jj);
		md[i]=(sx+sy)/2.0;
//		md[i]=sx;
//		if (sy<sx) md[i]=sy;		// pick the smallest length
	}

	// make a vector of arcs for new node location
//	vector_f nd;
	_ns.resize(0);
	float_t s=0;
	for (size_t i=0; i<n-1; i++)	// for each segment
	{
		float_t aver=(md[i]+md[i+1])/2.0;
		while (s<_curve.arcs().value_at(i+1))
		{
			push_back(_ns,s);			
			s+=aver;	
		}
	}
	if (_ns.size()<2) 
		throw "\nremesh fail in fjapp_bubbleNS.cpp\n";
	float_t ds=_curve.len()-_ns[_ns.size()-1];
	if (ds>0.125/npl)	// create extra one
	{
		s=(_ns[_ns.size()-2]+_curve.len())/2.0;
		_ns[_ns.size()-1]=s;
		push_back(_ns,_curve.len());
	}
	else
	_ns[_ns.size()-1]=_curve.len();
//	_curve.remesh(nd);
}

void TFJBubbleNS_Interface::calc_curvature(vector_f &cvt)
{
	size_t n=_curve.pt_count();
	cvt.resize(n);
	for (size_t i=0; i<n; i++)
		cvt[i]=(_xp2[i]*_yp1[i]-_yp2[i]*_xp1[i])/
					std::pow(_arc_crt[i],3.0);
	cvt[0]+=-_yp2[0]/_arc_crt[0];
	if (type_id==0)
		for (size_t i=1; i<n; i++)
			cvt[i]+=-_yp1[i]/_x[i]/_arc_crt[i];
	else {
		for (size_t i=1; i<n-1; i++)
			cvt[i]+=-_yp1[i]/_x[i]/_arc_crt[i];
		cvt[n-1]+=-_yp2[n-1]/_arc_crt[n-1];
	}

	TFJSmoothing sm_obj;
	sm_obj.set(16,16,4,true);
	vector_f tmp;
	for (size_t i=0; i<5; i++)
	{
		tmp=cvt;
		sm_obj.smoothing(tmp,cvt);
	}
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void TFJBubbleNS_Solver::initialize()
{
	// setup axis
	axis_type ax(2);

	float_t dx=1.0/_pms.xnpl, dy=1.0/_pms.ynpl;
/*
	float_t xpp0=0.1, dxn0=dx/4.0;
	float_t xpp1=0.5, dxn1=dx/2.0;
	ax.x().set_range(0,dxn0,xpp0/dxn0);
	ax.x().add_range(dxn1,(xpp1-xpp0)/dxn1);
	ax.x().add_range(dx,(_pms.xmax-xpp1)/dx);
	ax.y().set_range(0,dy,_pms.ymax/dy);
*/
//
	ax.x().set_range(0,dx,_pms.xmax/dx);
	ax.y().set_range(0,dy,_pms.ymax/dy);
//
	_axis.set(ax);

	// initialize surface
	initialize_interface();
//	_surf.set(_pms.snpl,_pms.flat);

	// initialize vof
	_vof.set_data(&_axis.axis(sptO),&_surf.curve());
	_vof.get_vof().eps=1e-12;
	_vof.set_Bo(_pms.Bo);
	_vof.set_coordinates(true);
#ifdef _Petsc_KSP
	initialize_ksp();
#endif
}

void TFJBubbleNS_Solver::reset_interface()
{ 
	if (_tc==0)
		_surf.reset(); 
	else
		spline_interface();		
}

void TFJBubbleNS_Solver::reset()
{
	// reset time
	_tc=_pms.t0;
	// reset surface 
//	if (_tc==0)
	reset_interface(); //_surf.reset();
//	else 
//		spline_interface(); // _surf.spline();
//	update_vof_st();
	// reset u,v,p
	_u.resize(_axis.xseg_count(sptU),_axis.yseg_count(sptU));
	_v.resize(_axis.xseg_count(sptV),_axis.yseg_count(sptV));
	_p.resize(_axis.xseg_count(sptP),_axis.yseg_count(sptP));
	if (_tc==0) {
		reset_matrix(_u); reset_matrix(_v); reset_matrix(_p);
	}
	// prepare for start
	solve_interface();
	prepare_step_scheme(false);
}

void TFJBubbleNS_Solver::step()
{
	// Solve NS eqns - 2nd order in time scheme
	solve_ns();
	solve_interface();
}

void TFJBubbleNS_Solver::after_step()
{
	prepare_step_scheme(true);
	_tc+=_pms.dt;
}

#ifndef _Petsc_KSP
void TFJBubbleNS_Solver::solve_ns()
{
	prepare_eqn_u();
	prepare_eqn_v();
//	_pp=_p;					// back up p for time n
	prepare_eqn_p();

	niter=0;
	bool Success=false,PSuccess=false,First=true;
	float_t acc_tmp=_pms.uvacc;
	do {								
//		if (stopped) return;

		Success|=PSuccess;		
		if (Success) acc_tmp=_pms.pacc;
		
		// STEP 1: Solve U hat

		// only update pressure part
		prepare_eqn_u(true);
		prepare_eqn_v(true);
		_pp=_p;
		CGSTAB(_axis.xseg_count(sptU),_axis.yseg_count(sptU),
				acc_tmp,_pms.max_iter2,_u,_uap,_uae,_uaw,_uan,_uas,_usu,_res);
		CGSTAB(_axis.xseg_count(sptV),_axis.yseg_count(sptV),
				acc_tmp,_pms.max_iter2,_v,_vap,_vae,_vaw,_van,_vas,_vsu,_res);

//		matrix_f _op;
		float_t _perr0=1.0,_perr;
		if (!Success)
		{
			// STEP 2: Calc U star, remove partial pressure gradient 
			correct_uv(1.0);
			// STEP 3: Solve P, new pressure gradient
//			_op=_p;		// OP is used for delta pressure check
			prepare_eqn_p(true);
			PSuccess=CGSTAB(_axis.xseg_count(sptP),_axis.yseg_count(sptP),
						_pms.pacc,_pms.max_iter2,_p,_pap,_pae,_paw,_pan,_pas,_psu,_res);
			// STEP 4: Calc U n+1, add partial pressure gradient back
			correct_uv(-1.0);

//			_op-=_p;
			_perr=mat_max_abs(_res);
			if (First) 
			{
				_perr0=_perr; _perr=1.0;
				First=false;
			}
			else
				_perr/=_perr0;
			PSuccess&=(_perr<_pms.acc);
		}
		if (PSuccess) Success=true;

		niter++;
	} while (!Success && (niter<_pms.max_iter1));	
}
#endif

void TFJBubbleNS_Solver::solve_interface()
{
	get_interface_uv();
	move_interface();
	prepare_interface();
}

void TFJBubbleNS_Solver::update_bc_u()
{
	int uxc=_axis.xseg_count(sptU);
	int uyc=_axis.yseg_count(sptU);
	for (int i=0; i<=uxc-1; i++) {
		_u(i,0)=_u(i,2)/3.0-2.0*_u(i,1);        // u=0 at bottom
		_u(i,uyc-1)=2*_u(i,uyc-2)-_u(i,uyc-3);	// d2u=0 at top
	}
	for (int j=0; j<=uyc-1; j++) {
		_u(0,j)=0;								// u=0 at center
		_u(uxc-1,j)=(4.0*_u(uxc-2,j)-_u(uxc-3,j))/3.0;	
												// 2nd order accuracy
	}
}

void TFJBubbleNS_Solver::update_bc_v()
{
//	float_t dx=_axis.dx(0);
//	float_t dy=_axis.dy(0);
//	int vxc=_axis.xseg_count();
//	int vyc=_axis.yseg_count();

	_axis.set_Stag(sptV);

	for (int i=1; i<=xc()-2; i++) {	
		// v=poiseuille flow at bottom
		if (xp(i)<_pms.IPos) 
			_v(i,0)=_pms.Uc*2.0*(1.0-xp(i)*xp(i));
		else _v(i,0)=0;			// v=0 at bottom
		
		int j=yc()-1;
		_v(i,j)=_v(i,j-1)-(_u(i,j)-_u(i-1,j))/dx(i)*dy(j)-
				(_u(i,j)+_u(i-1,j))/2.0/xp(i)*dy(j);
	}

	for (int j=0; j<yc(); j++) {
		int i=1;
		if (j!=0)									// coupled with u=0
			_v(i-1,j)=-2*_u(i,j)/dx(i)*dy(j)+_v(i,j-1);	// and continuity
		else
			_v(0,j)=_v(1,j);
		i=xc()-1;
		_v(i,j)=_v(i-1,j);					// dv/dx=0 at outter
	}
}

void TFJBubbleNS_Solver::correct_uv(float_t sign)
{
//	float_t dp,xp,xm;
//	float_t dx=_axis.dx(0);
//	float_t dy=_axis.dy(0);

//	int uxc=_axis.xseg_count(sptU);
//  int uyc=_axis.yseg_count(sptU);
	_axis.set_Stag(sptU);
    for (int i=1; i<xc()-1; i++) 	
		for (int j=1; j<yc()-1; j++) 
			_u(i,j)+=sign*_pms.Pn*_pms.beta*_pms.dt/_pms.Re/
					den(_uvof(i,j))/dx(i)*(_p(i+1,j)-_p(i,j));
	update_bc_u();

//    int vxc=_axis.xseg_count(sptV);
//    int vyc=_axis.yseg_count(sptV);
	_axis.set_Stag(sptV);
    for (int i=1; i<xc()-1; i++) 	
		for (int j=1; j<yc()-1; j++) 
			_v(i,j)+=sign*_pms.Pn*_pms.beta*_pms.dt/_pms.Re/
					den(_vvof(i,j))/dy(j)*(_p(i,j+1)-_p(i,j));
	update_bc_v();
}

void TFJBubbleNS_Solver::prepare_eqn_u(bool p_only)
{
//	float_t dx=_axis.dx(0);
//	float_t dx2=dx/2.0;
//	float_t dy=_axis.dy(0);

//	int uxc=_axis.xseg_count(sptU);
//  int uyc=_axis.yseg_count(sptU);

	_axis.set_Stag(sptU);

	// for each iteration, only need to change this part 
	if (p_only)
	{
	    for (int i=1; i<xc()-1; i++) 
			for (int j=1; j<yc()-1; j++) 
				_usu(i,j)+=_pms.Pn*(_pp(i+1,j)-_pp(i,j))/dx(i)-
							_pms.Pn*(_p(i+1,j)-_p(i,j))/dx(i);
		return;
	}

	_uap.resize(xc(),yc());	_usu.resize(xc(),yc());
	_uae.resize(xc(),yc());	_uaw.resize(xc(),yc());
	_uan.resize(xc(),yc());	_uas.resize(xc(),yc());
  
//	float_t xp,xe,xw,xm;
	float_t dp,ve,vw,vn,vs,vp;
//	float_t xtmp,ytmp;
//	float_t dxe,dxw,dxn,dxs;

//	ytmp=0.5/dy/dy;
    for (int i=1; i<xc()-1; i++) {
//		xp=_axis.axis(sptU).x().seg_center_value(i);
//		xe=xp(i)+dx2(i); 
//		xw=xp()-dx2(i);
//		dxe=(dx(i)+dx(i+1))/2.0;
//		dxw=(dx(i)+dx(i-1))/2.0;
//		xtmp=1.0/xp(i)/dx(i)/dx(i);
		for (int j=1; j<yc()-1; j++) 
		{
//			dye=(dx(i)+dx(i+1))/2.0;
//			dxw=(dx(i)+dx(i-1))/2.0;

			dp=den(_uvof(i,j));
			ve=vis(_pvof(i+1,j));	vw=vis(_pvof(i,j));
			vp=vis(_uvof(i,j));		
			vn=vis(_nvof(i,j));		vs=vis(_nvof(i,j-1));

			_uaw(i,j)=(xw(i)/xp(i))*vw/dx(i)/dxw(i);	//xtmp*xw*vw;
			_uae(i,j)=(xe(i)/xp(i))*ve/dx(i)/dxe(i);	//xtmp*xe*ve;
			_uan(i,j)=vn/2.0/dy(j)/dyn(j);	//ytmp*vn;
			_uas(i,j)=vs/2.0/dy(j)/dys(j);	//ytmp*vs;
			_uap(i,j)=-(_uaw(i,j)+_uae(i,j)+_uan(i,j)+_uas(i,j));

			_uap(i,j)=_uap(i,j)-dp*_pms.Re/_pms.dt-vp/xp(i)/xp(i);
					//-xtmp*(xe*ve+xw*vw)-vp/xp/xp-ytmp*(vn+vs); 

			_usu(i,j)=dp*_pms.Re*_pu(i,j)/_pms.dt+(3.0*_upcv(i,j)-_uppcv(i,j))/2.0-
					_pms.Pn*(_p(i+1,j)-_p(i,j))/dx(i)+_updf(i,j)/2.0+
					(3.0*_upst(i,j)-_uppst(i,j))/2.0/xp(i)/_pms.Ca/dx(i)/dy(j);
		}
    }
    for (int i=1; i<xc()-1; i++) 
	{
		_uap(i,1)+=-_uas(i,1);		_uas(i,1)=0;
		int j=yc()-2;
		_uap(i,j)+=_uan(i,j);		_uan(i,j)=0;
	}
	for (int j=1; j<yc()-1; j++)
	{
									_uaw(1,j)=0.0;
		int i=xc()-2;
		_uap(i,j)+=_uae(i,j);		_uae(i,j)=0;
	}
}

void TFJBubbleNS_Solver::prepare_eqn_v(bool p_only)
{
//	float_t dx=_axis.dx(0);
//	float_t dx2=dx/2.0;
//	float_t dy=_axis.dy(0);

//    int vxc=_axis.xseg_count(sptV);
//    int vyc=_axis.yseg_count(sptV);

	_axis.set_Stag(sptV);

	if (p_only)
	{
	    for (int i=1; i<xc()-1; i++) 
			for (int j=1; j<yc()-1; j++) 
				_vsu(i,j)+=_pms.Pn*(_pp(i,j+1)-_pp(i,j))/dy(j)-
							_pms.Pn*(_p(i,j+1)-_p(i,j))/dy(j);
		return;
	}

	_vap.resize(xc(),yc());	_vsu.resize(xc(),yc());
	_vae.resize(xc(),yc());	_vaw.resize(xc(),yc());
	_van.resize(xc(),yc());	_vas.resize(xc(),yc());
  
//	float_t xp(i),xe,xw,xmw,xme;
	float_t dp,ve,vw,vn,vs,vp;
//	float_t xtmp,ytmp;

//	ytmp=1.0/dy(j)/dy(j);
    for (int i=1; i<xc()-1; i++) {
//		xp(i)=_axis.axis(sptV).x().seg_center_value(i);
//		xe=xp(i)+dx(i)2; xw=xp(i)-dx(i)2;
//		xtmp=0.5/xp(i)/dx(i)/dx(i);
		for (int j=1; j<yc()-1; j++) {
			dp=den(_vvof(i,j));
			vn=vis(_pvof(i,j+1));	vs=vis(_pvof(i,j));
			ve=vis(_nvof(i,j));		vw=vis(_nvof(i-1,j));

			_vaw(i,j)=(xw(i)/xp(i))*vw/2.0/dx(i)/dxw(i);	//xtmp*xw*vw;
			_vae(i,j)=(xe(i)/xp(i))*ve/2.0/dx(i)/dxe(i);	//xtmp*xe*ve;
			_van(i,j)=vn/dy(j)/dyn(j);	//ytmp*vn;
			_vas(i,j)=vs/dy(j)/dys(j);	//ytmp*vs;
			_vap(i,j)=-(_vaw(i,j)+_vae(i,j)+_van(i,j)+_vas(i,j));

			_vap(i,j)=_vap(i,j)-dp*_pms.Re/_pms.dt;
				//-xtmp*(xe*ve+xw*vw)-ytmp*(vn+vs);

			_vsu(i,j)=dp*_pms.Re*_pv(i,j)/_pms.dt+(3.0*_vpcv(i,j)-_vppcv(i,j))/2.0-
					_pms.Pn*(_p(i,j+1)-_p(i,j))/dy(j)+_vpdf(i,j)/2.0+
					(3.0*_vpst(i,j)-_vppst(i,j))/2.0/xp(i)/_pms.Ca/dx(i)/dy(j);
	      }
    }
	float_t vv;
    for (int i=1; i<xc()-1; i++) 
	{
//		float_t pos=_axis.axis(sptV).x().seg_center_value(i);
		if (xp(i)<_pms.IPos) vv=_pms.Uc*2.0*(1.0-xp(i)*xp(i));
		else vv=0;
		_usu(i,1)+=_uas(i,1)*vv;		_uas(i,1)=0;
		int j=yc()-2;
		_uap(i,j)+=_uan(i,j);			_uan(i,j)=0;
	}
	for (int j=1; j<yc()-1; j++)
	{
		_uap(1,j)+=_uaw(1,j);			_uaw(1,j)=0.0;
		int i=xc()-2;
		_uap(i,j)+=_uae(i,j);			_uae(i,j)=0;
	}
}

void TFJBubbleNS_Solver::prepare_eqn_p(bool uv_only)
{
//	float_t dx=_axis.dx(0);
//	float_t dx2=dx/2.0;
//	float_t dy=_axis.dy(0);

//	int pxc=_axis.xseg_count(sptP);
//	int pyc=_axis.yseg_count(sptP);

	_axis.set_Stag(sptP);

	if (uv_only)
	{
//		float_t xp(i),xe,xw,
		float_t sutmp1,sutmp2;
	    for (int i=1; i<xc()-1; i++) 
		{
//			xp(i)=_axis.axis(sptP).x().seg_center_value(i);
//			xe=xp(i)+dx(i)2; xw=xp(i)-dx(i)2;
			for (int j=1; j<yc()-1; j++) 
			{
				sutmp1=sutmp2=0.0;
				if ((i!=1) && (i!=xc()-2))			// normal
					sutmp1=(xe(i)*_u(i,j)-xw(i)*_u(i-1,j))/dx(i)/xp(i);
				else
				{
					if (i==1)						// if at left
						sutmp1=xe(i)*_u(i,j)/dx(i)/xp(i);
				}
				if ((j!=1) && (j!=yc()-2))			// normal
					sutmp2=(_v(i,j)-_v(i,j-1))/dy(j);
				else
				{
			        if (j==1)						// if at bottom
					{
						sutmp2=_v(i,j)/dy(j);
						if (xp(i)<_pms.IPos) 
							sutmp2+=-_pms.Uc*2.0*(1.0-xp(i)*xp(i))/dy(j);
					}								// if at top
				}
				_psu(i,j)=(sutmp1+sutmp2)/_pms.Pn/_pms.beta;
			}
		}
		return;
	}

	_pap.resize(xc(),yc());	_psu.resize(xc(),yc());
	_pae.resize(xc(),yc());	_paw.resize(xc(),yc());
	_pan.resize(xc(),yc());	_pas.resize(xc(),yc());
  
//	float_t xp(i),xe,xw,xm;
//	float_t xtmpe,ytmpw;
	float_t ptmp1,ptmp2,sutmp1,sutmp2;
	float_t de,dw,dn,ds;
//	float_t uap0,uap1,vap0,vap1;

//	ytmp=1.0/dy(j)/dy(j);
	float_t reodt=_pms.Re/_pms.dt;
    for (int i=1; i<xc()-1; i++) 
	{
//		xp(i)=_axis.axis(sptP).x().seg_center_value(i);
//		xe=xp(i)+dx(i)2; xw=xp(i)-dx(i)2;
//		xtmp=1.0/xp(i)/dx(i)/dx(i); 
		for (int j=1; j<yc()-1; j++) 
		{
			_pap(i,j)=0;	_psu(i,j)=0;
			_paw(i,j)=0;	_pae(i,j)=0;
			_pas(i,j)=0;	_pan(i,j)=0;

			de=den(_uvof(i,j));	dw=den(_uvof(i-1,j));
			dn=den(_vvof(i,j));	ds=den(_vvof(i,j-1));

//			uap0=-de*reodt;	uap1=-dw*reodt;
//			vap0=-dn*reodt;	vap1=-ds*reodt;

			ptmp1=ptmp2=0.0;
			sutmp1=sutmp2=0.0;
			if ((i!=1) && (i!=xc()-2))			// normal
			{
				_pae(i,j)=-(xe(i)/xp(i))/dx(i)/dxe(i)/de;	//xtmp*xe/uap0;
				_paw(i,j)=-(xw(i)/xp(i))/dx(i)/dxw(i)/dw;	//xtmp*xw/uap1;
				ptmp1=-(_pae(i,j)+_paw(i,j));	//-xtmp*(xe/uap0+xw/uap1);
				sutmp1=(xe(i)*_u(i,j)-xw(i)*_u(i-1,j))/dx(i)/xp(i);
			}
			else
			{
				if (i==1)						// if at left
				{
					_pae(i,j)=-(xe(i)/xp(i))/dx(i)/dxe(i)/de;	//xtmp*xe/uap0;
					ptmp1=-_pae(i,j);	//-xtmp*xe/uap0;
					sutmp1=xe(i)*_u(i,j)/dx(i)/xp(i);
				}								// if at right
			}
			if ((j!=1) && (j!=yc()-2))			// normal
			{
				_pan(i,j)=-1.0/dy(j)/dyn(j)/dn;	//ytmp/vap0;
				_pas(i,j)=-1.0/dy(j)/dys(j)/ds;	//ytmp/vap1;
				ptmp2=-(_pan(i,j)+_pas(i,j));	//-ytmp*(1.0/vap0+1.0/vap1);
				sutmp2=(_v(i,j)-_v(i,j-1))/dy(j);
			}
			else
			{
                if (j==1)						// if at bottom
				{
					_pan(i,j)=-1.0/dy(j)/dyn(j)/dn;	//ytmp/vap0;
					ptmp2=-_pan(i,j);	//-ytmp/vap0;
					sutmp2=_v(i,j)/dy(j);
					if (xp(i)<_pms.IPos) 
						sutmp2+=-_pms.Uc*2.0*(1.0-xp(i)*xp(i))/dy(j);
				}								// if at top
			}
			_pae(i,j)/=reodt;	_paw(i,j)/=reodt;
			_pan(i,j)/=reodt;	_pas(i,j)/=reodt;
			_pap(i,j)=(ptmp1+ptmp2)/reodt;			
			_psu(i,j)=(sutmp1+sutmp2)/_pms.Pn/_pms.beta;

		}
    }
	_pap(xc()-2,yc()-2)=1.0;
}

void TFJBubbleNS_Solver::apply_interface_uv_bc()
{
	size_t cn=_surf.curve().pt_count();
	if (_surf.get_type()==0) {
		_us[0]=0;
		_us[cn-1]=0;	_vs[cn-1]=0;
	}
	else {
		_us[0]=0; 		_us[cn-1]=0;
	}
}

void TFJBubbleNS_Solver::get_interface_uv()
{
	update_bc_u();
	update_bc_v();

	size_t ii,jj;
	float_t nx,ny;

	size_t cn=_surf.curve().pt_count();
	_us.resize(cn);	_vs.resize(cn);
	_un.resize(cn);	_ut.resize(cn);
	for (size_t i=0; i<cn; i++)
	{
		float_t sx=_surf.curve().x(i);
		float_t sy=_surf.curve().y(i);

		// get interpolated u velocity
		ii=_axis.axis(sptV).x().locate_seg(sx);
		ii=bound(ii,(size_t)1,_axis.xseg_count(sptV)-2);
		nx=_axis.axis(sptV).x().seg_relative(ii,sx);
		jj=_axis.axis(sptV).y().locate_seg(sy);
		ny=_axis.axis(sptV).y().seg_relative(jj,sy);
		_us[i]=balance2(_u(ii-1,jj),_u(ii,jj),
						_u(ii,jj+1),_u(ii-1,jj+1),nx,ny);

		// get interpolated v velocity
		ii=_axis.axis(sptU).x().locate_seg(sx);
		nx=_axis.axis(sptU).x().seg_relative(ii,sx);
		jj=_axis.axis(sptU).y().locate_seg(sy);
		jj=bound(jj,(size_t)1,_axis.yseg_count(sptU)-2);
		ny=_axis.axis(sptU).y().seg_relative(jj,sy);
		_vs[i]=balance2(_v(ii,jj-1),_v(ii+1,jj-1),
						_v(ii+1,jj),_v(ii,jj),nx,ny);
	}

//	_us[0]=0;
//	_us[cn-1]=0;	_vs[cn-1]=0;
	apply_interface_uv_bc();

	// get normal and tangential velocity
	for (size_t i=0; i<cn; i++)
	{
		_un[i]=-_us[i]*_surf.yp(1)[i]+
					_vs[i]*_surf.xp(1)[i];
		_ut[i]=_us[i]*_surf.xp(1)[i]+
					_vs[i]*_surf.yp(1)[i];
	}
}

void TFJBubbleNS_Solver::move_interface()
{
	TFJBubbleNS_Interface::curve_type& cv=_surf.curve();
	for (size_t i=0; i<cv.pt_count(); i++)
	{
		cv.x(i)-=_un[i]*_surf.yp(1)[i]/
						_surf.arc_crt()[i]*_pms.dt;
		cv.y(i)+=_un[i]*_surf.xp(1)[i]/
						_surf.arc_crt()[i]*_pms.dt;
	}
}

void TFJBubbleNS_Solver::remesh_interface()
{
//	_surf.remesh(_axis.axis(sptO));
	_surf.remesh_uniform();
	_surf.finalize_mesh();
	_surf.spline();
}

void TFJBubbleNS_Solver::prepare_interface()
{
	spline_interface();
	remesh_interface();
/*
	for (size_t i=0; i<_surf.curve().pt_count(); i++)
		cout << i << " " << _surf.x()[i]
			<< "," << _surf.y()[i] << endl;
*/
//	_surf.x()[0]=-1e-5;
	update_vof_st();	
/*
	_surf.spline();
//	_surf.remesh_uniform();
	_surf.remesh(_axis.axis(sptO));
	_surf.finalize_mesh();
	_surf.spline();
//	if (_surf.remesh()) _surf.spline();
	update_vof_st();	
*/
}

void TFJBubbleNS_Solver::update_vof_st()
{
	_vof.vof_gen();	
	fill_vof_mat(sptU,_uvof);
	fill_vof_mat(sptV,_vvof);
	fill_vof_mat(sptP,_pvof);
	fill_vof_mat(sptN,_nvof);
	fill_st_mat(sptU,true,_ust);
	fill_st_mat(sptV,false,_vst);
}

bool TFJBubbleNS_Solver::CGSTAB(size_t xc, size_t yc,
							float_t acc, size_t niter,
							matrix_f& x, matrix_f& ap,
							matrix_f& ae, matrix_f& aw,
							matrix_f& an, matrix_f& as,
							matrix_f& b, matrix_f& res)
{
	res.resize(xc,yc);
	// Get initial residue of the eqn
	for (size_t i=1; i<xc-1; i++)
		for (size_t j=1; j<yc-1; j++)
		{
			res(i,j)=-ap(i,j)*x(i,j)-b(i,j)-
					ae(i,j)*x(i+1,j)-aw(i,j)*x(i-1,j)-
					an(i,j)*x(i,j+1)-as(i,j)*x(i,j-1);
		}
	matrix_f rhat=res;
	float_t err0=mat_max_abs(res);
	matrix_f pk(xc,yc),vk(xc,yc),sk(xc,yc),tk(xc,yc);
	reset_matrix(pk);	reset_matrix(vk);
	reset_matrix(sk);	reset_matrix(tk);
	float_t d=1,dp,alf=1,omeg=1,bet,tmp1,tmp2; 

	piter=0;
	// Start iteration
	float_t err=0;
	do
	{
		// calculate d
		dp=d;	d=0;
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				d+=rhat(i,j)*res(i,j);
		bet=d/dp*alf/omeg;

		// update p
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				pk(i,j)=res(i,j)+bet*(pk(i,j)-omeg*vk(i,j));

		// update v
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				vk(i,j)=ap(i,j)*pk(i,j)+
						ae(i,j)*pk(i+1,j)+aw(i,j)*pk(i-1,j)+
						an(i,j)*pk(i,j+1)+as(i,j)*pk(i,j-1);
		// alf
		tmp1=0;
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				tmp1+=rhat(i,j)*vk(i,j);
		alf=d/tmp1;

		// sk
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				sk(i,j)=res(i,j)-alf*vk(i,j);

		// tk
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				tk(i,j)=ap(i,j)*sk(i,j)+
						ae(i,j)*sk(i+1,j)+aw(i,j)*sk(i-1,j)+
						an(i,j)*sk(i,j+1)+as(i,j)*sk(i,j-1);

		// omeg
		tmp1=0;	tmp2=0;
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
			{
				tmp1+=tk(i,j)*sk(i,j);
				tmp2+=tk(i,j)*tk(i,j);
			}
		omeg=tmp1/tmp2;
		// x
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				x(i,j)+=alf*pk(i,j)+omeg*sk(i,j);		

		// res
		for (size_t i=1; i<xc-1; i++)
			for (size_t j=1; j<yc-1; j++)
				res(i,j)=sk(i,j)-omeg*tk(i,j);

//		err=mat_max_abs(res);
		err=(err0==0)?0:mat_max_abs(res)/err0;
		piter++;
	} while ((err>acc) && (piter<niter));	
	return (err<acc);
}

float_t TFJBubbleNS_Solver::mat_max_abs(matrix_f& ma)
{
	float_t mm=0.0;
	for (size_t i=1; i<ma.size1()-1; i++)
		for (size_t j=1; j<ma.size2()-1; j++)
		{
			float_t tmp=std::fabs(ma(i,j));
			if (tmp>mm) mm=tmp;
		}
	return mm;
}

void TFJBubbleNS_Solver::backup_data_u()
{
//	float_t dx=_axis.dx(0);
//	float_t dx2=dx/2.0;
//	float_t dy=_axis.dy(0);

	update_bc_u();
	_axis.set_Stag(sptU);
//    int uxc=_axis.xseg_
//    int uyc=_axis.yseg_count(sptU);
	_updf.resize(xc(),yc());
	_upst.resize(xc(),yc());	
	_upcv.resize(xc(),yc());	
  
//	float_t xp,xe,xw;
	float_t dp,ve,vw,vn,vs,vp;
//	float_t xtmp,ytmp;
//	ytmp=1.0/dy/dy;
    for (int i=1; i<xc()-1; i++) {
//		xp=_axis.axis(sptU).x().seg_center_value(i);
//		xe=xp+dx2; xw=xp-dx2;
//		xtmp=2.0/xp/dx/dx;
		for (int j=1; j<yc()-1; j++) {
			dp=den(_uvof(i,j));
			ve=vis(_pvof(i+1,j));	vw=vis(_pvof(i,j));
			vp=vis(_uvof(i,j));		
			vn=vis(_nvof(i,j));		vs=vis(_nvof(i,j-1));

			_updf(i,j)=((xe(i)/xp(i))*ve/dxe(i)*(_u(i+1,j)-_u(i,j))-
						(xw(i)/xp(i))*vw/dxw(i)*(_u(i,j)-_u(i-1,j)))*2.0/dx(i)
						-2.0*vp/xp(i)/xp(i)*_u(i,j)+
						(vn/dyn(j)*(_u(i,j+1)-_u(i,j))-
							vs/dys(j)*(_u(i,j)-_u(i,j-1)))/dy(j);

			_upcv(i,j)=-dp*_pms.Re*(_u(i,j)*(_u(i+1,j)-_u(i-1,j))/2.0/dx(i)+
				(_v(i,j)+_v(i+1,j)+_v(i,j-1)+_v(i+1,j-1))*
				(_u(i,j+1)-_u(i,j-1))/8.0/dy(j))+
				(vn*(_v(i+1,j)-_v(i,j))-vs*(_v(i+1,j-1)-_v(i,j-1)))/dx(i)/dy(j); ///2.0

			_upst(i,j)=_ust(i,j);
		}
    }
}

void TFJBubbleNS_Solver::backup_data_v()
{
//	float_t dx=_axis.dx(0);
//	float_t dx2=dx/2.0;
//	float_t dy=_axis.dy(0);

	update_bc_v();
	_axis.set_Stag(sptV);
//    int vxc=_axis.xseg_count(sptV);
//    int vyc=_axis.yseg_count(sptV);
	_vpdf.resize(xc(),yc());
	_vpst.resize(xc(),yc());	
	_vpcv.resize(xc(),yc());	

//	float_t xp(i),xe,xw;
	float_t dp,ve,vw,vn,vs,vp;
//	float_t xtmp,ytmp;

//	ytmp=2.0/dy(j)/dy(j);
    for (int i=1; i<xc()-1; i++) {
//		xp(i)=_axis.axis(sptV).x().seg_center_value(i);
//		xe=xp(i)+dx(i)2; xw=xp(i)-dx(i)2;
//		xtmp=1.0/xp(i)/dx(i)/dx(i);
		for (int j=1; j<yc()-1; j++) {
			dp=den(_vvof(i,j));
			vn=vis(_pvof(i,j+1));	vs=vis(_pvof(i,j));
			ve=vis(_nvof(i,j));		vw=vis(_nvof(i-1,j));

			_vpdf(i,j)=((xe(i)/xp(i))*ve/dxe(i)*(_v(i+1,j)-_v(i,j))
					-(xw(i)/xp(i))*vw/dxw(i)*(_v(i,j)-_v(i-1,j)))/dx(i)+
					(vn/dyn(j)*(_v(i,j+1)-_v(i,j))-
					vs/dys(j)*(_v(i,j)-_v(i,j-1)))*2.0/dy(j);

			_vpcv(i,j)=-dp*_pms.Re*((_u(i,j+1)+_u(i,j)+_u(i-1,j)+_u(i-1,j+1))
				*(_v(i+1,j)-_v(i-1,j))/8.0/dx(i)+_v(i,j)*
				(_v(i,j+1)-_v(i,j-1))/2.0/dy(j))+
				((xe(i)/xp(i))*ve*(_u(i,j+1)-_u(i,j))-
					(xw(i)/xp(i))*vw*(_u(i-1,j+1)-_u(i-1,j)))/dx(i)/dy(j); ///2.0;

			_vpst(i,j)=_vst(i,j);
	      }
    }
}

void TFJBubbleNS_Solver::prepare_step_scheme(bool second)
{
	if (!second)
	{ backup_data_u(); backup_data_v(); }
	_uppcv=_upcv;	_uppst=_upst;
	_vppcv=_vpcv;	_vppst=_vpst;
	if (second)
	{ backup_data_u(); backup_data_v(); }
	_pu=_u;	_pv=_v;	_pp=_p;
}

void TFJBubbleNS_Solver::fill_vof_mat(
						staggered_type stag_loc,
						matrix_f &mat)
{
	size_t m=_axis.xseg_count(stag_loc);
	size_t n=_axis.yseg_count(stag_loc);
	mat.resize(m,n);
	for (size_t i=0; i<m ;i++)
		for (size_t j=0; j<n; j++)
			mat(i,j)=_vof.get_staggered_vof(i,j,stag_loc);
}

void TFJBubbleNS_Solver::fill_st_mat(
						staggered_type stag_loc,
						bool fill_ust, matrix_f &mat)
{
	size_t m=_axis.xseg_count(stag_loc);
	size_t n=_axis.yseg_count(stag_loc);

	float_t u,v;
	mat.resize(m,n);
	for (size_t i=0; i<m ;i++)
		for (size_t j=0; j<n; j++)
		{
			_vof.get_staggered_st(i,j,stag_loc,u,v);
			if (fill_ust) mat(i,j)=u;
			else mat(i,j)=v;
		}
}

float_t TFJBubbleNS_Solver::get_tracked_vol()
{
	float_t _vol=0.0,ds;
	TFJPolynomial plx,ply,pl;
	for (size_t i=0; i<_surf.curve().seg_count(); i++)
	{
		ds=_surf.curve().arcs().seg_len(i);
		plx=_surf.curve().interp(0).get_seg_poly(i);
		ply=_surf.curve().interp(1).get_seg_poly(i);
		pl=(plx*plx*(ply>>1))<<1;
		_vol+=pl.eval_diff(0.0,1.0);
	}
	return -_vol*M_PI;
}

#ifdef _Petsc_KSP
void TFJBubbleNS_Solver::solve_ns()
{
	prepare_eqn_u();  	_uksp.initialize();
	prepare_eqn_v();	_vksp.initialize();
	_pp=_p;
	prepare_eqn_p();	 _pksp.initialize();

	// STEP 1: Solve U hat
	prepare_eqn_u(true); _uksp.solve();
	prepare_eqn_v(true); _vksp.solve();
	_pp=_p;

	// STEP 2: Calc U star, remove partial pressure gradient 
	correct_uv(1.0);
	// STEP 3: Solve P, new pressure gradient
	prepare_eqn_p(true); _pksp.solve();
	// STEP 4: Calc U n+1, add partial pressure gradient back
	correct_uv(-1.0);
	niter=1;
}

void TFJBubbleNS_Solver::initialize_ksp()
{
	_uksp.set_size(_axis.xseg_count(sptU),_axis.yseg_count(sptU));
	_vksp.set_size(_axis.xseg_count(sptV),_axis.yseg_count(sptV));
	_pksp.set_size(_axis.xseg_count(sptP),_axis.yseg_count(sptP));

	_uksp.set_data(&_uap,&_uaw,&_uae,&_uas,&_uan,&_usu,&_u);
	_vksp.set_data(&_vap,&_vaw,&_vae,&_vas,&_van,&_vsu,&_v);
	_pksp.set_data(&_pap,&_paw,&_pae,&_pas,&_pan,&_psu,&_p);

	_uksp.set_tolerance(_pms.uvacc);
	_vksp.set_tolerance(_pms.uvacc);
	_pksp.set_tolerance(_pms.pacc);
}
#endif

}		// end of namespace
