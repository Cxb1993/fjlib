
/*!
//	10.28 moved to seperat file
//			fixed dt counting error
//			add universal support for 
//			both surfactant and -free solver
//	12.19 added support for multiple nodes 
//	12.21 fixed a bug for dt_count due to integer division
//	12.30 fixed a bug for handling error "char const *"
//	12.31 add standardlized logging
//	2.24.06 changed to use fjapp_petsc.h, minimize petsc exposure
//			added support for diffusion, experiment
//	3.21.06 add support for rate_limit_step
!*/


#include "fjlib_log.h"
#include "fjlib_cio.h"
#include "fjapp_petsc.h"

// define a global varible log varible
#ifndef FJLIB_DISABLE_LOGGING
fjlib::TFJLog gLOG("log",&cout);
#endif

int	n_save,dt_count;

#include "fjlib_path.h"
#include "fjlib_vecmat_print.h"

using namespace fjlib;

TFJPetsc pet;

void save(my_Solver* obj)
{
	char *str=new char[256];
#ifdef __GNUC__ 
	gcvt(obj->tc(),8,str);
#else
	_gcvt(obj->tc(),8,str);
#endif
	str_t dn="t="+str_t(str);
	TFJPath::mkdir(dn);
	TFJPath::chdir(dn);

	quicksave(	"curve.txt",	obj->get_curve());
	save_mat(	"u.txt", 		obj->get_u_mat());
	save_mat(	"v.txt", 		obj->get_v_mat());
	save_mat(	"p.txt", 		obj->get_p_mat());
#ifdef _Surfactant
	save_vec(	"gama.txt",		obj->get_gama_vec());
#endif
#ifdef _Bulk_Diffusion
	save_mat(	"conc.txt",		obj->get_conc_mat());
	save_vec(	"fluxsub.txt",	obj->get_sublayer_flux());
	save_vec(	"cghost.txt",	obj->get_ghost_conc());
	save_vec(	"csub.txt",		obj->get_sublayer_conc());
#endif
	// for debug conveniency
	save_mat(	"vof.txt", 		obj->get_vof_mat());
	save_vec(	"xs.txt",		obj->get_surface().x());
	save_vec(	"ys.txt",		obj->get_surface().y());
	save_vec(	"su.txt",		obj->get_us_vec());	
	save_vec(	"sv.txt",		obj->get_vs_vec());
	save_vec(	"sun.txt",		obj->get_un_vec());
	save_vec(	"sut.txt",		obj->get_ut_vec());
	
	TFJPath::chdir("..");
}

void load(my_Solver* obj, const str_t& fn)
{
	TFJPath::chdir(fn);
	quickload("curve.txt", obj->get_curve());
	mprint_f ptr(&obj->get_u_mat());
	quickload("u.txt",ptr);
	ptr.set_matrix(&obj->get_v_mat());
	quickload("v.txt",ptr);
	ptr.set_matrix(&obj->get_p_mat());
	quickload("p.txt",ptr);
#ifdef _Surfactant
	vprint_f vptr(&obj->get_gama_vec());
	quickload("gama.txt", vptr);
#endif
#ifdef _Bulk_Diffusion
	ptr.set_matrix(&obj->get_conc_mat());
	quickload("conc.txt", ptr);
#endif

	TFJPath::chdir("..");	
}

#include "fjlib_string.h"
void update(void* obj)
{
	my_Solver* ob=(my_Solver*)obj;
	fjlib::float_t tc=ob->tc();

	gLOG << "t= " << tc << '(' << ob->get_iter1() << ')';
#ifndef _MULTIPLE_INTERFACE
	gLOG << "\tV: " << ob->get_tracked_vol() 
		<< "\tVe: " << ob->get_vol_err();
#else
	gLOG << "\tV: " << ob->get_Mvol_err();
#endif
	gLOG << cendl;
#ifdef _Surfactant
#ifndef _MULTIPLE_INTERFACE
	gLOG << "\t\tG: " << ob->get_curve_mass()
		<< "\tGe: " << ob->get_curve_mass_err();
#else
	gLOG << "\tM: " << ob->get_tracked_mass() << 
		'\t' << ob->get_total_Mmass();
#endif
#endif
	gLOG << cendl;
#ifdef _Bulk_Diffusion
	gLOG << "\t\tC: " << ob->get_bulk_mass()
		<< "\tCe: " << ob->get_bulk_mass_err();
#endif
	gLOG << cendl;

	dt_count++; 
	if (dt_count==n_save) { 
		save(ob); 
		dt_count=0;
		gLOG.SaveLog(); 
	}
}

#include "fjlib_memini.h"

#define IF_RANK_EQUAL_0 if(pet.rank==0)
#define IF_RANK_NOT_EQUAL_0 if(pet.rank!=0)

#include "mpi.h"

#define ROUND(x) (floor(x+0.5))

int main(int argc, char **args)
{
#ifdef _Petsc_KSP
	pet.Init(&argc,&args);
	IF_RANK_EQUAL_0
		gLOG << "using " << pet.size << " computer(s)" << cendl;
#endif
	double endtime,starttime;

	my_Solver solver;

	// make sure it doesn't save on other nodes
	IF_RANK_NOT_EQUAL_0
		gLOG.SetSaveOnDestroy(false);

	// NFS system for now
	my_Solver::params_type& pms=solver.get_params();
	TFJMemINI ini;

	quickload("bubble.ini",ini);
	pms.Re=ini.get_stream(1.0,"Re","constant");
	pms.Ca=ini.get_stream(0.3,"Ca","constant");
	pms.Bo=ini.get_stream(0.9,"Bo","constant");
	pms.Viscosity_Ratio=ini.get_stream(0.1,"viscosity_ratio","constant");
	pms.Density_Ratio=ini.get_stream(0.1,"density_ratio","constant");
	pms.Uc=ini.get_stream(1.0,"Uc","constant");
	pms.Pn=pms.IPos=1.0;
//	pms.Pn=pms.IPos=pms.Uc=1.0;

	pms.flat=ini.get_stream(true,"flat","geometry");
	pms.xmax=ini.get_stream(5.0,"xmax","geometry");
	pms.ymax=ini.get_stream(12.0,"ymax","geometry");
	pms.xnpl=ini.get_stream(12.0,"xnpl","geometry");
	pms.ynpl=ini.get_stream(12.0,"ynpl","geometry");
	pms.snpl=ini.get_stream(12.0,"snpl","geometry");

	pms.uvacc=ini.get_stream(1e-12,"uvacc","solver");
	pms.pacc=ini.get_stream(1e-12,"pacc","solver");
	pms.acc=ini.get_stream(1e-3,"acc","solver");
	pms.acc_abs=ini.get_stream(1e-2,"acc_abs","solver");
	pms.dt=ini.get_stream(1e-3,"dt","solver");
	// obsolete iters
	pms.max_iter1=ini.get_stream(2000,"outer_iters","solver");
	pms.max_iter2=ini.get_stream(500,"inner_iters","solver");

	fjlib::float_t dt_save=ini.get_stream(1e-1,"dt_save","main");
	n_save=ROUND(dt_save/pms.dt);
	pms.tmax=ini.get_stream(1.0,"tmax","main");
	pms.beta=1.0;

	pms.t0=ini.get_stream(0.0,"t0","main");
	str_t fn=ini.get_string("","t0_folder","main");
	if (pms.t0!=0)
		load(&solver,fn);
	gLOG << "solver info loaded." << cendl;
#ifdef _Surfactant
	my_Solver::params_surf_type& gpms=
							solver.get_surf_params();
	gpms.x0=ini.get_stream(0.1,"x","surfactant");
	gpms.Bi=ini.get_stream(0.1,"Bi","surfactant");
	gpms.Pe=ini.get_stream(1000,"Pe","surfactant");
	gpms.gacc=ini.get_stream(1e-4,"gacc","surfactant");
	gpms.Surfc.Vir=ini.get_stream(0.0,"Vir","surfactant");
	gpms.Surfc.El=ini.get_stream(0.2,"El","surfactant");
	gpms.Surfc.a=1.0;
	gpms.Surfc.GamaInf=1.0;
	gpms.Surfc.Sigma0=1.0;
	gLOG << "surfactant info loaded." << cendl;
#endif
#ifdef _Bulk_Diffusion
	my_Solver::params_diff_type& dpms=solver.get_diff_params();
	dpms.hob=ini.get_stream(1.0,"hob","diffusion");
	dpms.Peb=ini.get_stream(1.0,"Peb","diffusion");
	dpms.alpha=ini.get_stream(0.05,"alpha","diffusion");
	dpms.conc_acc=ini.get_stream(1e-12,"conc_acc","diffusion");
	dpms.iter_acc=ini.get_stream(1e-12,"iter_acc","diffusion");
	dpms.max_iters=ini.get_stream(1000,"max_iters","diffusion");
	int rls=ini.get_stream(3,"rate_limit","diffusion");
	dpms.rate_limit_step=TFJRateLimitStepType(rls);
//	bool diff=ini.get_stream(true,"switch_on","diffusion");
//	solver.switch_diffusion_on(diff);
#endif
#ifdef _MULTIPLE_INTERFACE
	gLOG << "multiple interface tracking enabled." << cendl;
#endif
	IF_RANK_EQUAL_0
 		solver.set(0,&update);
	solver.initialize();
	solver.reset();
//	return 0;

	IF_RANK_EQUAL_0 save(&solver);
	gLOG << "node " << pet.rank << " initialization done" << cendl;
	dt_count=0;
	IF_RANK_EQUAL_0 
		starttime=MPI::Wtime();
	try {
//		solver.step();
 		solver.run();
	}
	catch(char const *str) {		
		gLOG << "Caught exception: " << str << cendl;
		gLOG.SaveLog();
	}
	IF_RANK_EQUAL_0 {
		endtime=MPI::Wtime();
		gLOG << "Total time spent: " << endtime-starttime << "s" << cendl;
		gLOG << "Congratulations!! Done. " << cendl;
	}
}
