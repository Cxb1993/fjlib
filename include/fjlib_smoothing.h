#ifndef FJLIBSmoothingH
#define FJLIBSmoothingH

#include "fjlib_vecmat.h"

namespace fjlib {

/*!
//	Smoothing vector class
*/
class TFJSmoothing {
private:
	void		savgol(vector_f &c, int np, int nl, int nr,
						int ld, int m);
	void		ludcmp(matrix_f &a, vector_sz &indx, float_t &d);
	void		lubksb(matrix_f &a, vector_sz &indx, vector_f &b);
	float_t		householder_transform(vector_f &v);
	void		householder_hm(float_t tau, 
								const vector_f& v,
								matrix_f& a);
	void		qr_decomp(matrix_f& a,vector_f &tau);
	void		qr_unpack(const matrix_f& qr,
						const vector_f& tau,
						matrix_f& q, matrix_f& r);
	matrix_f	_q,_r;
	void		qr(const matrix_f& a,
					matrix_f& q, matrix_f& r,
					bool ecnomic=false);
protected:
    bool        with_ends;
    vector_f    sm_vec;
    int         n_left, n_right, n_order;
    void        calc_smvec();
	void		add_ends(const vector_f &c,
							vector_f &csm);
	vector_f	begin_pts(const vector_f &c);
	vector_f	end_pts(const vector_f &c);
public:
	TFJSmoothing(): with_ends(false) {}
	void		set(int left, int right, int order, bool ends=false);
    void		smoothing(const vector_f &c, vector_f& out);
	void		set_ends(bool _with_ends) 
	{ with_ends=_with_ends; }
};

}

#endif