#ifndef FJLIB_BubbleNS_MI_KSPH
#define FJLIB_BubbleNS_MI_KSPH

#include "fjapp_bubbleNS_MI.h"
#include "fjapp_ksp.h"

namespace fjlib {

class TFJBubbleNS_Solver_MI_KSP:public TFJBubbleNS_Solver_MI {
protected:
	TFJKSP			_uksp,_vksp,_pksp;
	void			solve_ns();
	int				_rank;
public:
	void			initialize();
	void			set_rank(int __rank) { _rank=__rank; }
	const
	int				rank() const { return _rank; }
};

}	// end of namespace

#endif
