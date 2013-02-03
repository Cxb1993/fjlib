//---------------------------------------------------------------------------
// Equation root finder wrapper

// Alas: interp

// To-do:
// 6.9.2005 revive it to new version with boost function
// 6.20		change the tempVec to be autosize

#ifndef FJLIBRootH
#define FJLIBRootH

#include "fjlib.h"
#include "fjlib_function.h"
#include <limits>

namespace fjlib {

// function form is not included in base due to the form uncentanty
class TFJRoot_Base {
protected:
	// Object where the data embeded
	void			*m_pObj;
	// Root valid range 
	float_t			m_fMin,m_fMax;
	// Solving tolerence
	float_t			m_fEps;
protected:
	vector_f		m_Roots;
	inline void		add_root(float_t x)
	{ push_back(m_Roots,x);	}
public:
	TFJRoot_Base(): m_pObj(NULL), 
					m_fEps(1e-12) {}
	inline
	size_t			root_count() { return m_Roots.size(); }
	virtual inline
	void			set_obj(void *obj)	{ m_pObj=obj; }
	void			set_range(float_t low, float_t high)
	{ m_fMin=low; m_fMax=high; }
	void			set_eps(float_t eps)  { m_fEps=eps; }

	virtual 
	void			solve() { m_Roots.resize(0); }
	bool			empty() { return (m_Roots.size()<1); }
	size_t			count() { return m_Roots.size(); }
	vector_f&		get_roots() { return m_Roots; }
};

class TFJRoot: public TFJRoot_Base {
protected:
	vector_f		m_TempVec;
	size_t			m_FuncOutParamNum;
	void			_auto_resize_vec(vector_f* v)
	{
		if (v->size()<m_FuncOutParamNum)
			v->resize(m_FuncOutParamNum);
	}
protected:
	TFJRoot(): m_FuncOutParamNum(1) 
	{
		_auto_resize_vec(&m_TempVec);
	}
	// Procedure which has 1 parameter in and N parameter out
	// Out vector, the first parameter is default storing function eval
	TFJProcIn1OutNp	m_Func;
	float_t			call(float_t x,vector_f* out)
	{ 
		_auto_resize_vec(out);
		m_Func(m_pObj,&x,out);
		return (*out)[0];
	}
	float_t			call_xf(float_t x)
	{
//		_auto_resize_vec(&m_TempVec);
		m_Func(m_pObj,&x,&m_TempVec);
		return m_TempVec[0];
	}
public:
	virtual
	void			set_func(TFJProcIn1OutNp func)
	{ m_Func=func; }

	float_t			root(int i) { return m_Roots[i]; }
};

class TFJRoot_Bracket: public TFJRoot {
	// Bracketing and Bisection method
protected:
	bool			brac(float_t& _x1,float_t& _x2);
	float_t			bisect(float_t _x1,float_t _x2);
public:
	void			solve()
	{ 
		TFJRoot::solve();
		float_t x1=m_fMin,x2=m_fMax;
		if (brac(x1,x2)) add_root(bisect(x1,x2)); 
	}
};

template<class T>
inline const T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

class TFJRoot_Brent: public TFJRoot {
	// Van Wijingaarden-Dekker-Brent Method
protected:
	float_t			brent();
public:
	void			solve()
	{ 
		TFJRoot::solve();
		add_root(brent()); 
	}
};

class TFJRoot_Newton: public TFJRoot {
	// Van Newton-Rapson Method using derivative 
protected:
	float_t			safe_newton(float_t x1, float_t x2);
	float_t			newton(float_t x1, float_t x2);
	bool			m_bSafe;
	void			call_xfdf(float_t x, float_t &f, float_t &df)
	{
		m_Func(m_pObj,&x,&m_TempVec);
		f=m_TempVec[0]; df=m_TempVec[1];
	}
public:
	TFJRoot_Newton(): m_bSafe(true) 
	{
		m_FuncOutParamNum=2;
		_auto_resize_vec(&m_TempVec);
	}
	void			set_safe(bool safe) { m_bSafe=safe; }
	void			solve()
	{ 
		TFJRoot::solve();
		if (m_bSafe) add_root(safe_newton(m_fMin,m_fMax));
		else add_root(newton(m_fMin,m_fMax)); 
	}
};

}	// end of namespace

#endif
