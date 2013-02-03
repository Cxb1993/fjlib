#ifndef FJAPP_BUBBLENSMI
#define FJAPP_BUBBLENSMI

#include "fjapp_bubbleNS.h"
#include "fjapp_MInterface.h"

namespace fjlib {

/*!
//
// The idea is to save work and use fjapp_bubbleNS 
// to handle interface as much as possible
// so whenever needed, data is exported from 
// TFJBubble_SInterface to local xy, when done,
// xy is imported into TFJBubble_SInterface for storage
//
// Cation:
// initialize_interface,reset_interface,solve_interface
// has been overwirrten for multiple interface handling 
// but other ???_interface rountine keep intact for 
// single interface handling !!! 
!*/
class TFJBubbleNS_Solver_MI: public TFJBubbleNS_Solver {
public:
	typedef TFJBubble_MInterface mi_type;
	typedef mi_type::iter_type mlist_iter_type;
protected:
	mi_type			_Msurf;			// MInterface object
	mlist_iter_type	_Ssurf;			// Current SInterface object

	vector_f		_Mus,_Mvs;		// Corresponding M storage objs
	vector_f		_Mun,_Mut;		
	vector_f		_Mcvt;		
	matrix_f		_Muvof,_Mvvof,_Mpvof,_Mnvof;
	matrix_f		_Must,_Mvst;

	vector_f		_Mvol;			// record vol 
private:
	// get rid of isolated node, the distance between 
	// whose two neighbours are almost zero
	void			_erase_mid_ridge();
	// get rid of end nodes which is already stick to center
	void			_erase_end_ridge();
	// get rid of small droplet
//	bool			_erase_droplet();
protected:
	virtual void	create_Minterface();

	// two recorded velocity used for repair Minterface
	float_t			_above_vel,_below_vel;	
	// check if split or merge is required and fix it
	void			repair_Minterface();

	virtual bool	split_Minterface();
	virtual bool	merge_Minterface();

	// refine and fetch a interface
	void			refine_Sinterface()
	{ _erase_mid_ridge(); _erase_end_ridge(); }
public:  // override this two major functions !!! hack
	// take care of multiple into pieces 
//	void			initialize_interface();
	void			reset_interface();
	void			solve_interface();
public:
	virtual void	import_interface(mlist_iter_type it);
	virtual void	export_interface(mlist_iter_type it);
public:
	float_t			get_tracked_Mvol() { return sum(_Mvol); }
	float_t			get_Mvol_err() 
	{ 
		float_t vol=get_total_vol();
		return (get_tracked_Mvol()-vol)/vol;
	}
};	// end of class

};	// end of namespace 

#endif