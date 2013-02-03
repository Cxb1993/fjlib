#include "fjapp_bubbleD.h"
#include "fjapp_petsc.h"
#include "fjlib_cio.h"
#include "fjlib_vecmat_print.h"
#include "fjlib_path.h"

using namespace fjlib;

TFJPetsc pet;
typedef TFJBubbleD_Solver my_Solver;

int n_save,dt_count;
#include "fjlib_string.h"
#include "fjlib_log.h"
TFJLog dlog("output.txt");

void save_data(my_Solver *obj)
{
	str_t dn="t="+to_string(obj->tc());
	TFJPath::mkdir(dn);
	TFJPath::chdir(dn);
	
	save_mat("map.txt",obj->get_unknown_map());
	quicksave("curve.txt",obj->get_curve());
	save_vec("csub.txt",obj->get_sublayer_conc());
	save_vec("cfluxsub.txt",obj->get_sublayer_cflux());
	save_mat("conc.txt",obj->get_bulk_conc());
	save_vec("cghost.txt",obj->get_ghost_conc());
	save_vec("su.txt",obj->get_su());
	save_vec("sv.txt",obj->get_sv());
	save_mat("vof.txt",obj->get_vof());
	save_mat("pmap.txt",obj->get_unknown_premap());

	TFJPath::chdir("..");
	
//	dlog << obj->tc() << " " << gavg << " " << cavg << cendl;
	dlog.SaveLog();

}

void update(void* obj)
{
	my_Solver* ob=(my_Solver*)obj;

	cout << "t=" << ob->tc() << "\t ";

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

#include "fjlib_memini.h"
int main(int argc, char *args[])
{
	pet.Init(&argc,&args);
	if (pet.rank==0)
		cout << "using " << pet.size << " computer(s)" << endl;
	
	my_Solver solver;
	my_Solver::params_type& p=solver.get_params();
	TFJMemINI ini;
	quickload("bubble.ini",ini);

	p.Pe=ini.get_stream(1.0,"Pe","constant");
	p.De=ini.get_stream(0.1,"De","constant");
	
	p.npl=ini.get_stream(11,"npl","geometry");
	p.snpl=ini.get_stream(11,"snpl","geometry");
	p.xmax=ini.get_stream(2.0,"xmax","geometry");
	p.ymax=ini.get_stream(5.0,"ymax","geometry");
	p.x0=ini.get_stream(1.0,"x0","geometry");
	p.ywell=ini.get_stream(0.0,"ywell","geometry");
	p.cylinder=ini.get_stream(0,"cylinder","geometry");
	
	p.dt=ini.get_stream(0.001,"dt","solver");
	p.t0=ini.get_stream(0.0,"t0","solver");
	p.tmax=ini.get_stream(1.0,"tmax","solver");
	p.conc_acc=ini.get_stream(1e-12,"acc","solver");
	n_save=ini.get_stream(100,"save_steps","solver");


/*
	p.Pe=pet.GetDouble("-Pe",1);
	p.npl=pet.GetInt("-npl",11);
	p.snpl=p.npl;
	p.xmax=pet.GetDouble("-xmax",2);
	p.ymax=pet.GetDouble("-ymax",5);
	p.dt=pet.GetDouble("-dt",0.001);
	p.t0=pet.GetDouble("-t0",0.0);
	p.tmax=pet.GetDouble("-tmax",1);
	p.conc_acc=1e-12;
	n_save=pet.GetInt("-nsave",100);
*/
	dt_count=0;
	if (pet.rank==0)
		solver.set(0,&update);

	if (pet.rank==0)
		cout << "initilizing ..." << endl;
	try {
		solver.initialize();
		solver.reset();
//		solver.step();
		solver.run();
	}
	catch (char const *str) {
		cout << "my error: " << str << endl;
	}
	if (pet.rank==0)
		cout << "Done." << endl;
}
