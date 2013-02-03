#include "fjapp_vof_staggered.h"
#include "fjlib_vecmat_print.h"
using namespace fjlib;

template <class ST,	class CT>
void fill_mat(TFJVOFGen_Staggered<ST,CT>& vof, 
			  TFJStaggeredPosType stag_loc,
				matrix_f &mat)
{
	TFJStaggeredAxis_Helper<ST> helper(*vof.get_axis());

	size_t m=helper.xseg_count(stag_loc);
	size_t n=helper.yseg_count(stag_loc);

	mat.resize(m,n);
	for (size_t i=0; i<m ;i++)
		for (size_t j=0; j<n; j++)
			mat(i,j)=vof.get_staggered_vof(i,j,stag_loc);
}

template <class ST,	class CT>
void fill_st(TFJVOFGen_Staggered_ST<ST,CT>& vof,				
			TFJStaggeredPosType stag_loc,
				matrix_f &mat)
{
	TFJStaggeredAxis_Helper<ST> helper(*vof.get_axis());
	size_t m=helper.xseg_count(stag_loc);
	size_t n=helper.yseg_count(stag_loc);

	fjlib::float_t u,v;
	mat.resize(m,n);
	for (size_t i=0; i<m ;i++)
		for (size_t j=0; j<n; j++)
		{
			vof.get_staggered_st(i,j,stag_loc,u,v);
			mat(i,j)=u;
		}
}

int main()
{
	typedef TFJCurve_CubicSpline curve_type;
	typedef TFJScale_Staggered_UniMesh scale_type;
	typedef TFJAxis_Staggered_UniMesh axis_type;
	typedef TFJVOFGen_Staggered_ST<scale_type,curve_type> vof_type;

	// set up axis
	float xnpl,ynpl,dx,dy,xmax,ymax;
	xnpl=ynpl=12.0;
	dx=1.0/xnpl; dy=1.0/ynpl;
	xmax=5.0; ymax=12.0;
	axis_type axis(2);
	axis.x().set_range(0,dx,xmax/dx);
	axis.y().set_range(0,dy,ymax/dy);

	// set up curve
	curve_type curve;
	vector_f x,y;
	curve.set_data(&x,&y);
	quickload("curve.txt",curve);
	curve.spline();

	// setup sigma
	vector_f sigma;
	vprint_f vptr(&sigma);
	quickload("sigma.txt",vptr);
	TFJInterp_Line sigInterp;
	sigInterp.set_data(curve.arcs().get_vector(),&sigma);
	
	// set up vof 
	vof_type vof;
	vof.set_data(&axis,&curve);
	vof.get_vof().eps=1e-12;
	vof.set_Bo(0.9);
	vof.set_coordinates(true);
	vof.set_STInterp(&sigInterp);

	// solve and fill matrix
	vof.vof_gen();
	matrix_f ust;
	fill_st(vof,sptU,ust);
	save_mat("ust.txt",ust);
}
