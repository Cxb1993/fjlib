#include "fjapp_diffGCC.h"
#include "fjlib_cio.h"
#include "fjlib_vecmat_print.h"

namespace fjlib {

void TFJDiffGC_Coupled::initialize_interp()
{
	_gc.set_env(&_genv);
	_gc.gen_env_info(_gc.mirror_pos(),2,1);
}

void TFJDiffGC_Coupled::prepare_interp()
{
	vector_n src(3);
	src[0]=0; src[1]=1; src[2]=2;
	_gc.fill_env_nodes(_conc,_csub);
	_gc.gen_interp_coef(src,_gcoef);
}

void TFJDiffGC_Coupled::interp_gnodes()
{
	size_t ii,jj,x,y;

	float_t tmp;
	// for each ghost, update values
	gc_type::npos_vec_type& nv=_gc.ghosts_npos();
	gc_type::pos_vec_type& mp=_gc.mirror_pos();
	vector_n& mnp=_gc.msurf_npos();
	for (size_t i=0; i<nv.size(); i++)
	{
		ii=nv[i].x();	jj=nv[i].y();
		tmp=_gc.interp_xy(i,mp[i].x(),mp[i].y(),_gcoef);
		_conc(ii,jj)=2.0*_csub[mnp[i]]-tmp;	// assume c=1 on surface
//		_conc(ii,jj)=_gc.interp_xy(i,x,y,_gcoef);
	}
}

void TFJDiffGC_Coupled::initialize()
{
	TFJDiffGC::initialize();
	_csub.resize(_surf.pt_count());
	reset_vector(_csub,1.0);
}

float_t TFJDiffGC_Coupled::_check_err(const matrix_f& p,
									const matrix_f& q)
{
	float_t err=0,tmp;
	for (size_t i=0; i<p.size1(); i++)
		for (size_t j=0; j<p.size2(); j++)
		{
			tmp=fabs(p(i,j)-q(i,j));
			if (tmp>err) { err=tmp; _ii=i; _jj=j; }
		}
//	print_var("ii",ii);
//	print_var("jj",jj);
	return err;
}

void TFJDiffGC_Coupled::solve_Dirichlet()
{
//	print_mat("umat",_gc.get_umat(),vmpfText | vmpfColMajor | vmpfColReverse); 
	// doing intialization for the first time
	initialize_interp();
/*	
	for (size_t i=0; i<_genv.size(); i++)
	{
		gc_type::node_vec_type& nv=_genv[i];
		cout << i << " ";
		for (size_t j=0; j<nv.size(); j++)
			cout << "( " << nv[j].i << "," << nv[j].j << ") ";
		cout << endl;
	}
*/
	prepare_interp();
	interp_gnodes();
//	print_mat("coef",_gcoef,vmpfText); 
	prepare();
	_cksp.initialize();
	_cksp.solve();		// solve
//	print_mat("conc",_conc,vmpfText | vmpfColMajor | vmpfColReverse); 
	// from now, it can be solved couple of times
	matrix_f cp;
	float_t err=10;
	_iters=1;
	while (err>_tol) {
		cp=_conc;
		prepare_interp();
		interp_gnodes();	// update ghost nodes
		prepare(true);		// prepare again
		_cksp.solve();
//		print_mat("conc",_conc,vmpfText | vmpfColMajor | vmpfColReverse); 
		err=_check_err(cp,_conc);
		print_var("err",err);
		_iters++;
	}
}

void TFJDiffGC_Coupled::solve()
{
	solve_Dirichlet();
}

}	// end of namespace
