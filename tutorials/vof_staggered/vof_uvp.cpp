#include "fjlib_test.h"

#include "fjlib_float.h"
using namespace fjlib;

static float_t MY_FLOAT_EPS=FJLIB_FLOAT_EPS*1000;

#include "fjapp_vof_staggered.h"
#include "fjlib_vecmat_print.h"
void print_mat(matrix_f& m)
{
	mprint_f mp(&m,vmpfText|vmpfColMajor|vmpfColReverse);
	std::cout << mp << std::endl;
}

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

struct testDataFixture1 {
	vector_f x,y;
	TFJCurve_Line cv;
	testDataFixture1() 	{
		set_vector(x,0.0,0.25,0.75,0.75,1.5,1.0,0.0);
		set_vector(y,1.5,0.75,0.75,1.5 ,1.5,0.5,0.0);
		cv.set_data(&x,&y);
	}
};

struct testUniMeshFixture {
	testDataFixture1 d;
	TFJAxis_Staggered_UniMesh m;
	TFJVOFGen_Staggered<TFJScale_Staggered_UniMesh,
					TFJCurve_Line> gen;

	testUniMeshFixture() { 
		m.redim(2);
		m.x().set_range(0,0.5,2);
		m.y().set_range(0,0.5,2);
		gen.set_data(&m,&d.cv);
		gen.get_vof().eps=MY_FLOAT_EPS;
	}
};

struct testMeshFixture {
	testDataFixture1 d;
	TFJAxis_Staggered_Mesh m;
	TFJVOFGen_Staggered<TFJScale_Staggered_Mesh,
					TFJCurve_Line> gen;

	testMeshFixture() { 
		m.redim(2);
		m.x().set_range(0,0.5,2);
		m.y().set_range(0,0.5,2);
		gen.set_data(&m,&d.cv);
		gen.get_vof().eps=MY_FLOAT_EPS;
	}
};

struct testMeshFixture2 {
	testDataFixture1 d;
	TFJAxis_Staggered_Mesh m;
	TFJVOFGen_Staggered<TFJScale_Staggered_Mesh,
					TFJCurve_Line> gen;

	testMeshFixture2() { 
		m.redim(2);
		vector_f x;
		set_vector(x,0.0,0.25,1.0);
		m.x().set_values(x);
		m.y().set_range(0,0.5,2);
		gen.set_data(&m,&d.cv);
		gen.get_vof().eps=MY_FLOAT_EPS;
	}
};

struct testSTFixture {
	testDataFixture1 d;
	TFJAxis_Staggered_Mesh m;
	TFJVOFGen_Staggered_ST<TFJScale_Staggered_Mesh,
					TFJCurve_Line> gen;

	testSTFixture() { 
		m.redim(2);
		m.x().set_range(0,0.5,2);
		m.y().set_range(0,0.5,2);
		gen.set_data(&m,&d.cv);
		gen.get_vof().eps=MY_FLOAT_EPS;
	}
};

void test_UniMesh()
{
	testUniMeshFixture tf;
	tf.gen.vof_gen();
	matrix_f vm;
	fill_mat(tf.gen,sptO,vm);

//	print_mat(tf.gen.get_vofmat());
	print_mat(vm);

	BOOST_CHECK_CLOSE(vm(0,1),17.0/24.0,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(1,1),0.75,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(1,0),0.25,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(0,0),0.75,MY_FLOAT_EPS);
}

void test_Mesh()
{
	testMeshFixture tf;
	tf.gen.vof_gen();
	matrix_f vm;
	fill_mat(tf.gen,sptO,vm);
	print_mat(vm);
//	print_mat(tf.gen.get_vofmat());

	BOOST_CHECK_CLOSE(vm(0,1),17.0/24.0,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(1,1),0.75,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(1,0),0.25,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(0,0),0.75,MY_FLOAT_EPS);

/*
	std::cout << std::endl << std::endl;
	vector_f r;
	r=tf.gen.get_vof().roots_at(3,3,bptS);
	std::cout << r << std::endl;
	r=tf.gen.get_vof().roots_at(3,3,bptE);
	std::cout << r << std::endl;
	r=tf.gen.get_vof().roots_at(3,3,bptN);
	std::cout << r << std::endl;
	r=tf.gen.get_vof().roots_at(3,3,bptW);
	std::cout << r << std::endl;
*/
}

void test_Mesh2()
{
	testMeshFixture2 tf;
	tf.gen.vof_gen();
	matrix_f vm;
	fill_mat(tf.gen,sptO,vm);
	print_mat(vm);
//	print_mat(tf.gen.get_vofmat());

	BOOST_CHECK_CLOSE(vm(0,1),11.0/12.0,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(1,1),2.0/3.0,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(1,0),3.0/8.0,MY_FLOAT_EPS);
	BOOST_CHECK_CLOSE(vm(0,0),7.0/8.0,MY_FLOAT_EPS);
}

void test_ST()
{
	testSTFixture tf;
	tf.gen.vof_gen();
	matrix_f vm;
	fill_st(tf.gen,sptO,vm);
	print_mat(vm);
}

test_suite*	fjlib_init_test()
{
    test_suite* test= BOOST_TEST_SUITE( "FJScale Test" );	
    test->add( BOOST_TEST_CASE( &test_dummy));
    test->add( BOOST_TEST_CASE( &test_UniMesh));
    test->add( BOOST_TEST_CASE( &test_Mesh));
    test->add( BOOST_TEST_CASE( &test_Mesh2));
    test->add( BOOST_TEST_CASE( &test_ST));

    return test;
}

BOOST_TEST_MAIN_MACRO


/*
void fill_mat(TFJVOFGen_Staggered_ST& vof, 
			  TFJStaggeredPosType stag_loc,
						matrix_f &mat)
{
	TFJStaggered_UniMesh2D mesh(*vof.get_axis());
	size_t m=mesh.xseg_count(stag_loc);
	size_t n=mesh.yseg_count(stag_loc);

	mat.resize(m,n);
	for (size_t i=0; i<m ;i++)
		for (size_t j=0; j<n; j++)
			mat(i,j)=vof.get_staggered_vof(i,j,stag_loc);
}

void fill_st(TFJVOFGen_Staggered_ST& vof,				
			TFJStaggeredPosType stag_loc,
						matrix_f &mat)
{
	TFJStaggered_UniMesh2D mesh(*vof.get_axis());
	size_t m=mesh.xseg_count(stag_loc);
	size_t n=mesh.yseg_count(stag_loc);

	float_t u,v;
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
	TFJVOFGen::axis_type axis(2);
	TFJVOFGen::curve_type curve;

	axis.x().set_range(0,0.5,2);
	axis.y().set_range(0,0.5,2);

	vector_f x(7),y(7);
	x[0]=0;			y[0]=1.5;
	x[1]=0.25;		y[1]=0.75;
	x[2]=0.75;		y[2]=0.75;
	x[3]=0.75;		y[3]=1.5;
	x[4]=1.5;		y[4]=1.5;
	x[5]=1;			y[5]=0.8;
	x[6]=0;			y[6]=0;
	curve.set_data(&x,&y);
	curve.spline();
	cout << "Length= " << curve.len() << endl;
//	cout << curve.interp(1).get_seg_poly(5) << endl;

	TFJVOFGen_Staggered_ST vf;
	vf.set_data(&axis,&curve);
	vf.vof_gen();

	mprint_f mf(&vf.get_vofmat(),vmpfText|vmpfColMajor|vmpfColReverse);
	cout << mf << endl;
	quicksave("refined.txt",mf);
//	test_cell(vf,0,0);

//	vf.set_coordinates(true);
	matrix_f mat;
	mprint_f smf(&mat,vmpfText|vmpfColMajor|vmpfColReverse);

	cout << "OOO" << endl;
	fill_mat(vf,sptO,mat);
	cout << smf << endl;

	cout << "PPP" << endl;
	fill_mat(vf,sptP,mat);
	cout << smf << endl;
	cout << "UUU" << endl;
	fill_mat(vf,sptU,mat);
	cout << smf << endl;
	cout << "VVV" << endl;
	fill_mat(vf,sptV,mat);
	cout << smf << endl;
	cout << "NNN" << endl;
	fill_mat(vf,sptN,mat);
	cout << smf << endl;

	quicksave("new.txt",smf);

	// testing surface tension
	matrix_f ust;
	fill_st(vf,sptP,ust);
	mprint_f mfst(&ust,vmpfText|vmpfColMajor|vmpfColReverse);
	cout << mfst << endl;


	//  tesing the original mesh comparison
	TFJVOFGen_Full vff;
	vff.set_data(&axis,&curve);
	vff.vof_gen();
	matrix_f matf;
	vff.calc_vofs(matf);
	vff.fill_vofs(matf);
	mprint_f smff(&matf,vmpfText|vmpfColMajor|vmpfColReverse);
	cout << smff << endl;
	quicksave("old.txt",smff);

	int j;
	cin >> j;
}
*/
