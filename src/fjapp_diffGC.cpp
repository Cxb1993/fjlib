#include "fjapp_diffGC.h"
//#include "fjlib_vecmat_print.h"

namespace fjlib {

void TFJDiffGC::initialize()
{
	// create a geometry coordinates, half spherical x>0
	float_t dlen=1.0/_npl;;
	int nn=int(_nlen*_npl+0.5);
	_axis.redim(2);
	_axis.x().set_range(0.0,dlen,nn);
	_axis.y().set_range(-_nlen,dlen,nn*2);
//	cout << "x range " << _axis.x().min() << "-" << _axis.x().max() << endl;
//	cout << "y range " << _axis.y().min() << "-" << _axis.y().max() << endl;
	
	// create a half spherical surface
	size_t nc=M_PI*_bradius*_npl+1;
	if (nc<5) throw "come on, give more nodes in fjapp_diffGC.cpp";
	_x.resize(nc); _y.resize(nc);
	float_t dsita=M_PI/(nc-1);
	float_t ac;
	for (size_t i=0; i<nc; i++)
	{
		ac=dsita*i;
		_x[i]=_bradius*std::sin(ac);
		_y[i]=_bradius*std::cos(ac);
	}
	_x[0]=-1e-5;
	_x[nc-1]=-1e-5;
//	cout << "x coord: " << endl << _x << endl;
//	cout << "y coord: " << endl << _y << endl;

	//initial matrix and vector size
	size_t m=_axis.x().seg_count(), n=_axis.y().seg_count();
//	cout << "matrix size: " << m << "x" << n << endl;
	_cap.resize(m,n);	_csu.resize(m,n);
	_caw.resize(m,n);	_cae.resize(m,n);
	_cas.resize(m,n);	_can.resize(m,n);
	_conc.resize(m,n);	_pre_umat.resize(m,n);

	// initialize surf
	_surf.set_data(&_x,&_y);
	typedef TFJCubicSpline at;
	_surf.interps().x().set_bc(at::csbtClamped,
							at::csbtClamped,1,-1);
	_surf.interps().y().set_bc(at::csbtClamped,
							at::csbtClamped,0,0);
	_surf.calc_arcs();
	_surf.spline();

	// setup vof
	_vof.set_data(&_axis,&_surf);
	_vof.set_coordinates(true);
	_vof.vof_gen();
	_vof.fill_unknowns(_pre_umat);
//	print_mat("pre_umat",_pre_umat,vmpfText | vmpfColMajor | vmpfColReverse);
	
	// setup gc
	_gc.set_data(&_axis,&_surf,_pre_umat);
	_gc.gc_gen();
//	print_mat("umat",_gc.get_umat(),vmpfText | vmpfColMajor | vmpfColReverse);

	// initial parallel
#ifdef _Petsc_KSP
	initialize_ksp();
#endif
}

void TFJDiffGC::_set_ghosts()
{
	// for each ghost, update values
	gc_type::npos_vec_type& nv=_gc.ghosts_npos();
	for (size_t i=0; i<nv.size(); i++)
		_conc(nv[i].x(),nv[i].y())=1.0;
}

void TFJDiffGC::initialize_ksp()
{
	_cksp.set_data(&_cap,&_caw,&_cae,
					&_cas,&_can,&_csu,&_conc);
	_cksp.set_mask(&_gc.get_umat());
	_cksp.set_tolerance(_acc);
}

void TFJDiffGC::prepare(bool ghost_only)
{
	if (!ghost_only)
	{
	size_t m=_pre_umat.size1(), n=_pre_umat.size2();
	
	// generate coefficient for bulk
	reset_matrix(_cap,-4.0);
	reset_matrix(_caw,1.0);
	reset_matrix(_cae,1.0);
	reset_matrix(_cas,1.0);
	reset_matrix(_can,1.0);
	reset_matrix(_csu,0.0);

	// correction for cylindrical coordinates
	float_t dlen=_axis.x().seg_len(0);	// assume uniform mesh
	float rp,rt;
	int nn=int(_nlen*_npl+0.5);
	for (size_t i=0; i<m; i++)
	{
		rp=_axis.x().seg_center_value(i);
		rt=dlen/2.0/rp;
		for (size_t j=0; j<n; j++)
		{
			_cae(i,j)=1+rt;
			_caw(i,j)=1-rt;
		}	
	}
	
	// generate coefficient for regular b.c
	size_t i,j;
	// left and right
	for (j=0; j<n; j++)
	{
		i=0;	_caw(i,j)=_cas(i,j)=_can(i,j)=0.0;	_cae(i,j)=4.0;
//		i=m-1;	_cae(i,j)=_cas(i,j)=_can(i,j)=0.0;	_caw(i,j)=4.0;
		i=m-1;	_cae(i,j)=_cas(i,j)=_can(i,j)=_caw(i,j)=_csu(i,j)=0.0;
	}
	// top and bottom
	for (i=1; i<m-1; i++)
	{
		j=0;	_cae(i,j)=_cas(i,j)=_can(i,j)=_caw(i,j)=_csu(i,j)=0.0;
		j=n-1;	_cae(i,j)=_cas(i,j)=_can(i,j)=_caw(i,j)=_csu(i,j)=0.0;
//		j=0;	_cas(i,j)=_caw(i,j)=_cae(i,j)=0.0;	_can(i,j)=4.0;
//		j=n-1;	_can(i,j)=_caw(i,j)=_cae(i,j)=0.0;	_cas(i,j)=4.0;
	}
	// set top right corner
//	i=0; j=0;	_cae(i,j)=0.0;
	
	}	// end of ghost_only
	
	// generate coefficient for ghost b.c
	gc_type::npos_vec_type& pv=_gc.ghosts_npos();
	size_t ii,jj;
	size_t k=0;
	while (k<pv.size()) {
		ii=pv[k].x(); jj=pv[k].y();
		_caw(ii,jj)=_cae(ii,jj)=_cas(ii,jj)=_can(ii,jj)=0.0;
		_csu(ii,jj)=4.0*_conc(ii,jj);
		k++;
	}
/*	
	print_mat("caw",_caw,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("cae",_cae,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("cas",_cas,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("can",_can,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat("csu",_csu,vmpfText | vmpfColMajor | vmpfColReverse);
*/
}

void TFJDiffGC::solve()
{
	_set_ghosts();
	prepare();
	_cksp.initialize();
	_cksp.solve();
}

}	// end of namespace
