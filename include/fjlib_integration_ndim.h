#ifndef FJINTEGRATEM_H
#define FJINTEGRATEM_H

// 1.4 integrate with fixed dim value is added

#include "fjlib_integration.h"

namespace fjlib {

/*!
//	Multi-dimension Integration class
//
//	Function form is no different from TFJIntegrate, therefore
//	value at different dimensions has to be accessed through
//	data object TFJIntegrate_Data and its derivative
//
//	Call integrate_fix() for dim-locked integration, 
//	ex. integration f(x,y) from x=0 to x=2, with y fixed at 3
//
//
//	1.4	integrate with certain fixed dim is added
*/
class TFJIntegrate_NDim: public TFJIntegrate_Base {
public:
	typedef TFJIntegrate* integ_ptr;
	typedef std::vector<integ_ptr> integ_vec;
	///	Total dimension of the integration
	int			dim;
	///	Vector of 1D integration objects
	integ_vec
				integs;
	///	Set dimension
	void		set_dimension(int d) { dim=d; integs.resize(d); }
	/// Default constructor, 1D
	TFJIntegrate_NDim() { set_dimension(1); }
	/// Constructor, nD
	TFJIntegrate_NDim(int d) { set_dimension(d); }
public:
	///	Default setup call
	void		set(TFJIntegrate_Data* data_obj,
						func_type func, bool reset=true)
	{
		TFJIntegrate_Base::set(data_obj,func,reset);
		setup();
	}
protected:
	virtual
	void		setup()
	{
		TFJIntegrate_NDim b;
		func_type func=boost::bind(
					&TFJIntegrate_NDim::func_ndim,&b,_1,_2);
		fixed.resize(dim);
		fixed_value.resize(dim);
		for (size_t i=0; i<dim; i++)
			fixed[i]=false;
		if (dim>1)
			for (int i=0; i<dim-1; i++)
				integs[i]->set(data,func,i+1);
		// user defined function
		integs[dim-1]->set(data,function,dim);
	}
public:
	float_t		value_at(float_t* x)
	{ 
		for (size_t i=0; i<dim; i++)
			data->x[i]=x[i];
		return function(data,x);
	}
	float_t		integrate()
	{ 		
		dim_index=1;
		return integs[0]->integrate(); 
	}
protected:
	/// Info of fixed dimension
	vector_b	fixed;
	/// Info of fixed value at each dimension
	vector_f	fixed_value;
public:
	/// Integrate fixed dim call
	virtual
	float_t		integrate_fix(const vector_b& f,
								const vector_f& v)
	{
		fixed=f; fixed_value=v;
		dim_index=0;
		while ((dim_index<dim) && f[dim_index])
		{
			data->x[dim_index]=v[dim_index];
			dim_index++;
		}
		dim_index++;
		if (dim_index>dim)
			return value_at(&(data->x[0]));
		else
			return integs[dim_index-1]->integrate();
	}
public:
	float_t		func_ndim(void* obj, float_t* x)
	{
		float_t result;
		TFJIntegrate_Data *f=(TFJIntegrate_Data *)obj;
		TFJIntegrate_NDim *integ=
			(TFJIntegrate_NDim *)f->integs;
		int i=integ->dim_index;
		f->x[i-1]=*x;

		while ((integ->dim_index<integ->dim) &&
				integ->fixed[integ->dim_index])
		{
			f->x[integ->dim_index]=
				integ->fixed_value[integ->dim_index];
			integ->dim_index++;
		}
		integ->dim_index++;
		if (integ->dim_index>integ->dim)
			result=integ->value_at(&(f->x[0]));
		else
			result=integ->integs[integ->dim_index-1]->integrate();

		integ->dim_index=i;
		return result;
	}
};


} // end of namespace

#endif