#include "fjlib_cio.h"
#include "fjapp_gc.h"
#include "fjlib_vecmat_print.h"

using namespace fjlib;

struct testDataFixture1 {
	vector_f x,y;
	TFJCurve_Line cv;
	TFJAxis_UniMesh me;
	testDataFixture1() 	{
		set_vector(x,0.0,0.50,0.75,1.0,0.75);
		set_vector(y,1.5,1.25,1.00,0.5,0.00);
//		set_vector(x,0.0,0.25,0.75,0.75,1.5,1.0,0.0);
//		set_vector(y,1.5,0.75,0.75,1.5 ,1.5,0.5,0.0);
		cv.set_data(&x,&y); 

		me.redim(2);
		me.x().set_range(0,0.25,8);
		me.y().set_range(0,0.25,8);   
}
};

   
template <class ST,class CT>
void setup_vof_gc(TFJVOFGen_Full<ST,CT>& vf, matrix_f& mat, matrix_n& mat2)
{
//	vf.eps=1e-7;
	vf.vof_gen();
	vf.calc_vofs(mat);
	vf.fill_vofs(mat);
	vf.fill_unknowns(mat2);
}
   
struct testGeoFixture1 {
//private:
	typedef TFJGCGen<TFJScale_UniMesh,TFJCurve_Line> gc_type;
	testDataFixture1 df;
	TFJVOFGen_Full<TFJScale_UniMesh,TFJCurve_Line> vof;
	matrix_f vofm;
	matrix_n uknm;
public:
	gc_type gc;

	testGeoFixture1() { 
		vof.set_data(&df.me,&df.cv);
		setup_vof_gc(vof,vofm,uknm);
		gc.set_data(&df.me,&df.cv,uknm);
	}
}; 

void testGC1() {
	testGeoFixture1 gf;
	print_mat(gf.vofm,vmpfText | vmpfColMajor | vmpfColReverse);
	print_mat(gf.uknm,vmpfText | vmpfColMajor | vmpfColReverse);
 
	gf.gc.gc_gen();     
	print_mat(gf.gc.get_umat(),vmpfText | vmpfColMajor | vmpfColReverse);
	
	print_var("ghosts count",gf.gc.ghosts_count());
	print_var("bnodes count",gf.gc.bnodes_count());

	for (size_t i=0; i<gf.gc.ghosts_count(); i++)
		cout << "(" << gf.gc.ghosts_npos()[i].x() << ","
				<<	gf.gc.ghosts_npos()[i].y() << ")" << endl;

	typedef testGeoFixture1::gc_type	 gc_type;

	// test search algoritm here
	double x,y;
	x=gf.df.me.x().seg_center_value(gf.gc.ghosts_npos()[0].x());
	y=gf.df.me.y().seg_center_value(gf.gc.ghosts_npos()[0].y());
	gc_type::node_vec_type nv;
	gf.gc.search_grid_nodes(x,y,nv,3);
	for (size_t i=0; i<nv.size(); i++)
		cout << nv[i].x << "," << nv[i].y 
			<< "[" << nv[i].i << "," << nv[i].j << "]D"
			<< nv[i].dist << endl;
	cout << endl;

	// test search mirror node here
	for (size_t i=0; i<gf.gc.ghosts_count(); i++)
	{
		gc_type::pos_type &mp=gf.gc.mirrors_pos()[i];
//		gc_type::pos_type &mgp=gf.gc.mghosts_pos()[i];
		cout << i << ":";
		cout << " s=" << gf.gc.mirrors_arc()[i] 
			<< " m(" << mp.x() << "," << mp.y() << ")" << endl;
//			g(" 
//			<< mgp.x() << "," << mgp.y() << ")" << endl;
	}
	

	// make a ghosts nodes list
	
/*
	TFJGCGen<TFJScale_UniMesh,TFJCurve_Line>::node_vec_type nv;
	size_t n=7;
	float x=0.625, y=1.25;
	gf.gc.search_grid_nodes(x,y,nv,n);
	for (size_t i=0; i<n; i++)
		cout << nv[i].x << "," << nv[i].y 
			<< "[" << nv[i].i << "," << nv[i].j << "]D"
			<< nv[i].dist << endl;
	cout << endl;
	nv.resize(0);
	n=3;
	gf.gc.search_surf_nodes(x,y,nv,n);
	for (size_t i=0; i<n; i++)
		cout << nv[i].x << "," << nv[i].y 
			<< "[" << nv[i].i << "," << nv[i].j << "]D"
			<< nv[i].dist << endl;
*/
//	gf.gc.gen_ghosts_env();
} 

#include "fjlib_test.h"
 
test_suite*	fjlib_init_test()
{
    test_suite* test= BOOST_TEST_SUITE( "GC Test" );	
    test->add( BOOST_TEST_CASE( &test_dummy));
	test->add( BOOST_TEST_CASE( &testGC1));
    return test;
}
  	  
BOOST_TEST_MAIN_MACRO
/*
int main()
{ 
	cout << "come on" << endl;
	testGC1();	
	cout << "" << endl;  
}
*/
