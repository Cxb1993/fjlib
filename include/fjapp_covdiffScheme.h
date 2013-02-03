#ifndef FJAPP_Convection_Diffusion_Scheme
#define FJAPP_Convection_Diffusion_Scheme

#include "fjlib.h"

namespace fjlib {

/*! 
// only convection term is evaluated here
// diffusion term remains same
!*/
class TFJCovDiffScheme {
protected:
	vector_f	*_f;
	vector_f	*_v;
	virtual
	void		process();
	// interpolated value for positive and negtive part
	vector_f	fpos,fneg;
private:
	virtual
	float_t		HatD2N(float_t val,
						float_t upstream,
						float_t downstream);
	virtual
	float_t		HatN2D(float_t val,
						float_t upstream,
						float_t downstream);
protected:
	virtual 
	float_t		HatH2F(float_t v) { return v; }
public:
	void		set_data(vector_f* val, vector_f* vel)
	{ _f=val; _v=vel; process(); }
	void		get_west(size_t i, float_t& cp, 
						float_t& cw, float_t& cr);
	void		get_east(size_t i, float_t& cp,
						float_t& ce, float_t& cr);
};

class TFJCovDiffScheme_MINMOD: public TFJCovDiffScheme {
protected:
	float_t		HatH2F(float_t v);
};

class TFJCovDiffScheme_SMART: public TFJCovDiffScheme {
protected:
	float_t		HatH2F(float_t v);
};

class TFJCovDiffScheme_STOIC: public TFJCovDiffScheme {
protected:
	float_t		HatH2F(float_t v);
};

class TFJCovDiffScheme_HOAB: public TFJCovDiffScheme {
protected:
	float_t		HatH2F(float_t v);
};

class TFJCovDiffScheme_USED: public TFJCovDiffScheme {
protected:
	float_t		HatH2F(float_t v);
};

}	// end of namespace 

#endif