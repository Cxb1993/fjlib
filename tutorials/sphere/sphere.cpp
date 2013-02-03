#include "fjlib_cio.h"
#include "fjapp_SphereAdsorp.h"
#include "fjapp_petsc.h"
#include "fjlib_vecmat_print.h"
#include "fjlib_path.h"
 
using namespace fjlib;

TFJPetsc pet;
typedef TFJSphereAdsorp my_Solver;
int n_save,dt_count;

#include "fjlib_string.h"
#include "fjlib_log.h"
TFJLog dlog("output.txt");

void save_data(my_Solver *obj)
{
	str_t dn="t="+to_string(obj->tc());
	TFJPath::mkdir(dn);
	TFJPath::chdir(dn);
	
	save_vec("gama.txt",obj->gama());
	save_mat("map.txt",obj->unknown_map());
	save_vec("csub.txt",obj->conc_sublayer());
	save_vec("cfluxsub.txt",obj->cflux_sublayer());
	save_mat("conc.txt",obj->conc_bulk());

	TFJPath::chdir("..");
	
	// save average of gama and csub
	size_t n=obj->gama().size();
	float gavg,cavg;
	gavg=norm_1(obj->gama())/n;
	cavg=norm_1(obj->conc_sublayer())/n;
	dlog << obj->tc() << " " << gavg << " " << cavg << cendl;
	dlog.SaveLog();

}

void update(void* obj)
{
	my_Solver* ob=(my_Solver*)obj;

	cout << "t=" << ob->tc() << "\t[ ";
	for (size_t i=0; i<ob->oiter_count(); i++)
		cout << ob->iiter_counts()[i] << " ";
	cout << "]";

	dt_count++;
	if (dt_count==n_save) 
	{
		save_data(ob);
		dt_count=0;
		cout << " saved." << endl;
	}
	else
		cout << endl;
}

int main(int argc, char *args[])
{
	pet.Init(&argc,&args);
	if (pet.rank==0)
		cout << "using " << pet.size << " computer(s)" << endl;
	
	my_Solver sphere;
	my_Solver::params_type& p=sphere.get_params();

	p.hob=pet.GetDouble("-hob",20);
	p.Dsb=pet.GetDouble("-Dsb",0.01);
	p.x_eq=pet.GetDouble("-xeq",0.697);
	p.lamda=pet.GetDouble("-lamda",10.71);
//	p.x0=pet.GetDouble("-x0",0.0);
	p.radius=pet.GetInt("-r",1);
	p.npl=pet.GetInt("-npl",2);
	p.size=pet.GetInt("-size",5);
	p.dt=pet.GetDouble("-dt",0.001);
	p.alpha=pet.GetDouble("-alpha",0.05);
	p.conc_acc=1e-12;
	p.iter_acc=pet.GetDouble("-acc",1e-12);
	p.debug=pet.GetInt("-debug",0);
	p.max_iters=500;
	p.tmax=pet.GetDouble("-tmax",p.dt);
	n_save=pet.GetInt("-nsave",1);
	dt_count=0;
	if (pet.rank==0)
		sphere.set(0,&update);

	if (pet.rank==0)
		cout << "initilizing ..." << endl;
	try {
		sphere.initialize();
		sphere.reset();
//		sphere.step();
		sphere.run();
	}
	catch (char const *str) {
		cout << "my error: " << str << endl;
	}
//	matrix_f& s=sphere.conc_bulk();
/*
	if (s.size1()<12) 
		print_mat("solution",s,vmpfText | vmpfColMajor | vmpfColReverse);
	else {
		save_mat("out.txt",s);
		cout << "solution saved to out.txt" << endl;
	}
*/
	if (pet.rank==0)
		cout << "Done." << endl;
}
