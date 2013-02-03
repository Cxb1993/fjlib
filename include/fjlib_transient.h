#ifndef TFJLIBTransientH
#define TFJLIBTransientH

#include "fjlib.h"
#include "fjlib_function.h"

namespace fjlib {


/*!
//	Framework for transient process
// created 2002 probably
// 2.21.2006 add after_step() func
*/
class TFJTransient {
protected:
	/// Callback before step function
	TFJProcedure	beforestep_func;
	/// Callback after step function
	TFJProcedure	afterstep_func;
	void			on_afterstep()
	{ if (!afterstep_func.empty()) afterstep_func(this); }
	void			on_beforestep()
	{ if (!beforestep_func.empty()) beforestep_func(this); }
	/// Internal storage for steps count
	size_t			step_count;
	/// Internal trigger for termination
	bool			stopped;
public:
	void			set(const TFJProcedure& before, 
						const TFJProcedure& after)
	{ beforestep_func=before; afterstep_func=after; }
	///	Return step count
	inline
	size_t			count() { return step_count; }
	/// Setup one time intialization work before use
	virtual
	void			initialize() 
	{ step_count=0; }
	/// Reset to the initial condition
	virtual
	void			reset()
	{ on_afterstep(); }
	///	Run till stopped preset or manually
	virtual
	void			run()
	{
		// example
		stopped=false;
		while (!stopped && !is_end()) {
			on_beforestep();		// can be inside step() 
			step();
			after_step();
			on_afterstep();
			step_count++;		
		};
	}
	///	Step to the next
	virtual
	void			step()=0;
	virtual
	void			after_step() {}
	/// Stop the process manually
	virtual
	void			stop() 
	{ stopped=true; }
	///	Condition for terminate the process 
	virtual
	bool			is_end()=0;
};

}

#endif
