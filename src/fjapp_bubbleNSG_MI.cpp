#include "fjapp_bubbleNSG_MI.h"

// 1.10.06 changed according to fjapp_bubbleNS_MI.cpp
//		thinking about merging this into fjapp_bubbleNS.cpp
//		but realize that fjapp_bubbleNSG are derived from
//		fjapp_bubbleNS and they are too different.
//		therefore have to mantain two almost identical
//		fjapp_bubbleNS?_MI.cpp files, whihc is lame
//		One possible solution is to derive everything from
//		fjapp_bubbleNS_MI.cpp to fjapp_bubbleNS_MIG.cpp

#include "fjlib_log.h"
extern fjlib::TFJLog gLOG;

namespace fjlib {


///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
void TFJBubbleNSG_Solver_MI::create_Minterface()
{
	_Msurf.clear();
	// copy local to SInterface
	mlist_iter_type it=_Msurf.append(),it1;
	export_interface(it);

	// split them into pieces if needed
	bool found=true;
	while (found)
	{
		matrix_f& mat=it->data();
		found=false;
		for (size_t i=0; i<mat.size2(); i++)
		{
			if (mat(drX,i)<0) {
				found=true;
				// split at this location
				it1=_Msurf.append();
				it->split(*it1,i+1,false);
				it=it1;				
			}
		}	
	}
}

void TFJBubbleNSG_Solver_MI::reset_interface()
{
	// if t=0, reset interface
	if (_tc==0)
	{
		_surf.reset();
		size_t cn=_surf.curve().pt_count();
		_gama.resize(cn);
		reset_vector(_gama,1.0);
	}
	// no matter what, create Minterface
	create_Minterface();
	_tracked_mass=-1;
}

void TFJBubbleNSG_Solver_MI::import_interface(mlist_iter_type it)
{
	it->export_pos(_surf.x(),_surf.y());
	it->export_row(drGama,_gama);
	mlist_iter_type e=it; e++;
	_surf.set_type((e!=_Msurf.mlist().end()));
	spline_interface();	
}

void TFJBubbleNSG_Solver_MI::export_interface(mlist_iter_type it)
{
	it->import_pos(_surf.x(),_surf.y());
	it->import_row(drGama,_gama);
}

void TFJBubbleNSG_Solver_MI::solve_interface()
{
//	size_t cn=_Msurf.pt_size(); 
	_Mus.resize(0);	_Mvs.resize(0);
	_Mun.resize(0);	_Mut.resize(0);

	// initial value
	_above_vel=-100;
	_Mvol.resize(0);
	_Mmass.resize(0);
	
	// hack for tracked mass, first time
	bool ft=false;
	if (_tracked_mass<0) { ft=true; _tracked_mass=0; }

	// iterate each interface
	mlist_iter_type cb=_Msurf.mlist().begin();
	for (_Ssurf=cb; _Ssurf!=_Msurf.mlist().end(); _Ssurf++) 
	{
		// repair ridge if needed
		refine_Sinterface();
		// first dump to local surf
		import_interface(_Ssurf);
		// do normal stuff 
		TFJBubbleNSG_Solver::solve_interface();
		// fill data needed back to the global varaibles
		export_interface(_Ssurf);
		push_back(_Mvol,get_tracked_vol());
		push_back(_Mmass,get_total_mass());
		_tracked_mass+=get_absorbed_mass();

		// 3. add up all the vof st
		if (_Ssurf==cb) {
			_Muvof=_uvof; _Mvvof=_vvof;
			_Mpvof=_pvof; _Mnvof=_nvof;
			_Must=_ust;   _Mvst=_vst;
		} else {
			_Muvof+=_uvof; _Mvvof+=_vvof;
			_Mpvof+=_pvof; _Mnvof+=_nvof;
			_Must+=_ust;   _Mvst+=_vst;
		}
		// handle split merge of interface
		_below_vel=_un[0];
		repair_Minterface();
		_above_vel=_un[_un.size()-1];
	}
	_uvof=_Muvof; _vvof=_Mvvof;
	_pvof=_Mpvof; _nvof=_Mnvof;
	_ust=_Must; _vst=_Mvst;

	if (ft) _tracked_mass+=get_total_Mmass();

	// export Minterface for output, hack
	_Msurf.export_pos(_surf.x(),_surf.y());	
	_Msurf.export_row(drGama,_gama);
}

void TFJBubbleNSG_Solver_MI::repair_Minterface()
{
	bool handled=split_Minterface();
	if (!handled)
	{
		float_t rel_vel=_below_vel+_above_vel;
		if (rel_vel>0)
			handled=merge_Minterface();
	}
}

bool TFJBubbleNSG_Solver_MI::split_Minterface()
{
	// check if split is needed
	bool split=false;
	int i;
	int cn=_surf.curve().pt_count();
	if ((_Msurf.size()<10) && (cn>8))
	{
		int mind=5;
		i=mind;
		while (i<cn-mind)
		{
			if (_surf.x()[i]<1.0/_surf.get_npl()) {
				split=true; break;
			}
			i++;
		}
	}
	if (split)	
	{
		// print info about splitting
		gLOG << "split i=" << i << cendl;
		gLOG << "split after x=" << _surf.x()[i] << " y=" << _surf.y()[i] << cendl;
		gLOG << "next one x=" << _surf.x()[i+1] << " y=" << _surf.y()[i+1] << cendl;
		// insert one after _Ssurf
		mlist_iter_type itt=_Ssurf; itt++;
		mlist_iter_type it=_Msurf.insert(itt);
		_Ssurf->split(*it,i+1,false);
		// switch to the new surf
		_Ssurf=it;
	}
	return split;
}

bool TFJBubbleNSG_Solver_MI::merge_Minterface()
{
	bool merge=false;
	// no merge if it's the starting interface
	if (_Ssurf==_Msurf.mlist().begin()) return merge;

	// find y of both interface
	mlist_iter_type it=_Ssurf; it--;
	matrix_f& m0=it->data(), m1=_Ssurf->data();
	float_t above_y=m0(drY,m0.size2()-1),
			below_y=m1(drY,0);

	if (fabs(above_y-below_y)<1.0/_surf.get_npl())
	{
		gLOG << "merge between y=" << above_y << " and " << below_y << cendl;
		it->merge(*_Ssurf,1);
		_Msurf.mlist().erase(_Ssurf);
		// back to the previous surf
		_Ssurf=it;
		merge=true;
	}
	return merge;
}

void TFJBubbleNSG_Solver_MI::_erase_mid_ridge()
{
	float_t dx,dy,ds;
	matrix_f &mat=_Ssurf->data();

	int i=1;
	bool found=false;
	while (!found && (i<_Ssurf->size()-1))
	{
		dx=mat(drX,i+1)-mat(drX,i-1); 
		dy=mat(drY,i+1)-mat(drY,i-1);
		ds=std::sqrt(dx*dx+dy*dy);
		if (ds<0.5/_surf.get_npl()) found=true;
		else i++;
	}
	if (found)	// remove one node from the list
		_Ssurf->erase(i,i);
}

void TFJBubbleNSG_Solver_MI::_erase_end_ridge()
{
	using namespace std;
	matrix_f &mat=_Ssurf->data();

	// check the beginning of the curve
	int i=0;
	while (mat(drX,i)<1e-2)
	{
		float_t dx=mat(drX,i+1)-mat(drX,i);
		if (fabs(dx)>1e-2) break;
		i++;
	}
	if (i!=0)	// take the latter part
		_Ssurf->erase(0,i);

	i=_Ssurf->size()-1;
	while (mat(drX,i)<1e-2)
	{
		float_t dx=mat(drX,i-1)-mat(drX,i);
		if (fabs(dx)>1e-2) break;
		i--;
	}
	if (i!=_Ssurf->size()-1)	// take the first part
		_Ssurf->resize(i+1);
}

}	// end of namespace
