#include "fjlib_cio.h"
#include "fjlib_io.h"
#include "fjlib_vecmat_print.h"
#include "qplot_matrix.h"

int main()
{
	qplot::TQPMatrix a;
	fjlib::quickload("in.txt",a);
	fjlib::print_mat(a.mat);
}
