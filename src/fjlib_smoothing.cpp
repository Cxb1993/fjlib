#include "fjlib_smoothing.h"

namespace fjlib {

void TFJSmoothing::savgol(vector_f &c, int np, int nl, int nr,
							int ld, int m)
{
	int j,k,imj,ipj,kk,mm;
	float_t d,fac,sum;

	if (np < nl+nr+1 || ld > m || nl+nr < m)
		throw "\nerror in savgol() in fjlib_smoothing.cpp\n";
	vector_sz indx(m+1);
	matrix_f a(m+1,m+1);
	vector_f b(m+1);
	for (ipj=0;ipj<=(m << 1);ipj++) {
		sum=(ipj ? 0.0 : 1.0);
		for (k=1;k<=nr;k++) sum += std::pow(float_t(k),float_t(ipj));
		for (k=1;k<=nl;k++) sum += std::pow(float_t(-k),float_t(ipj));
		mm=std::min(ipj,2*m-ipj);
		for (imj = -mm;imj<=mm;imj+=2) a[(ipj+imj)/2][(ipj-imj)/2]=sum;
	}
	ludcmp(a,indx,d);
	for (j=0;j<m+1;j++) b[j]=0.0;
	b[ld]=1.0;
	lubksb(a,indx,b);
	for (kk=0;kk<np;kk++) c[kk]=0.0;
	for (k = -nl;k<=nr;k++) {
		sum=b[0];
		fac=1.0;
		for (mm=1;mm<=m;mm++) sum += b[mm]*(fac *= k);
		kk=(np-k) % np;
		c[kk]=sum;
	}
}

void TFJSmoothing::ludcmp(matrix_f &a, vector_sz &indx, float_t &d)
{
	const float_t TINY=1.0e-20;
	int i,imax,j,k;
	float_t big,dum,sum,temp;

	int n=a.size1();
	vector_f vv(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=std::fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) 
			throw "\nsingular matrix in routine ludcmp in fjlib_smoothing.cpp\n";
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*std::fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}

void TFJSmoothing::lubksb(matrix_f &a, vector_sz &indx, vector_f &b)
{
	int i,ii=0,ip,j;
	float_t sum;

	int n=a.size1();
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum != 0.0)
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void TFJSmoothing::set(int left, int right, int order, bool ends)
{
	n_left=left;
	n_right=right;
	n_order=order;
	with_ends=ends;
	calc_smvec();
}

void TFJSmoothing::calc_smvec()
{
    int n=n_left+n_right+1;
    sm_vec.resize(n);
	savgol(sm_vec,n,n_left,n_right,0,n_order);
}

void TFJSmoothing::smoothing(const vector_f &c, vector_f& out)
{
	size_t cn=c.size();
	if (cn<n_order) { out=c; return; }
    int n=n_left+n_right+1;
	if (cn<n)
	{
		n=cn;					// switch to less span
		div_t dr=div(n,2);
		if (dr.rem==0) n=n-1;	// make it a odd number
		n_left=(n-1)/2;
		n_right=n_left;
		calc_smvec();
	}
    out=c;
	for (size_t i=n_left; i<cn-n_right; i++)
	{
		// plug in
		float_t tmp=0.0;
		for (int j=-n_left; j<=n_right; j++)
        {
			if (j<0) tmp+=sm_vec[n+j]*c[i+j];
			else	tmp+=sm_vec[j]*c[i+j];
        }
		out[i]=tmp;
	}
	if (with_ends) add_ends(c,out);
}

float_t TFJSmoothing::householder_transform(vector_f &v)
{
	using namespace std;
	size_t n=v.size();
	if (n<=1) return 0.0;
	float_t xnorm=0.0;
	for (size_t i=1; i<n; i++)
		xnorm+=v[i]*v[i];
	xnorm=sqrt(xnorm);
	if (xnorm==0.0) return 0.0;
	float_t alpha=v[0],
			hypot=sqrt(alpha*alpha+xnorm*xnorm),
			beta,tau;
	if (alpha>=0) beta=-hypot;
	else beta=hypot;
	tau=(beta-alpha)/beta;
	for (size_t i=1; i<n; i++)
		v[i]/=alpha-beta;
	v[0]=beta;
	return tau;
}

void TFJSmoothing::householder_hm(float_t tau, 
							const vector_f& v,
							matrix_f& a)
{
	if (tau==0) return;
	size_t m=a.size1(),
			n=a.size2();
	float_t wj;
	for (size_t j=0; j<n; j++)
	{
		wj=a(0,j);
		for (size_t i=1; i<m; i++)
			wj+=a(i,j)*v[i];
		a(0,j)-=tau*wj;
		for (size_t i=1; i<m; i++)
			a(i,j)-=tau*v[i]*wj;
	}
}

void TFJSmoothing::qr_decomp(matrix_f& a,vector_f &tau)
{
	size_t m=a.size1(),
			n=a.size2();
	size_t d=std::min(m,n);
	tau.resize(d);
	float_t tau_i;
	for (size_t i=0; i<d; i++)
	{
		vector_f v(m-i);
		for (size_t j=i; j<m; j++)
			v[j-i]=a(j,i);
		tau_i=householder_transform(v);
		for (size_t j=i; j<m; j++)
			a(j,i)=v[j-i];

		tau[i]=tau_i;
		if (i+1<n)
		{
			matrix_f mm(m-i,n-i-1);
			for (size_t ii=i; ii<m; ii++)
				for (size_t jj=i+1; jj<n; jj++)
					mm(ii-i,jj-i-1)=a(ii,jj);
			householder_hm(tau_i,v,mm);
			for (size_t ii=i; ii<m; ii++)
				for (size_t jj=i+1; jj<n; jj++)
					a(ii,jj)=mm(ii-i,jj-i-1);
		}
	}
}

void TFJSmoothing::qr_unpack(const matrix_f& qr,
						const vector_f& tau,
						matrix_f& q, matrix_f& r)
{
	size_t m=qr.size1(),
			n=qr.size2();
	size_t d=std::min(m,n);
	q.resize(m,m);
	r.resize(m,n);
	for (size_t i=0; i<m; i++)
        for (size_t j=0; j<m; j++)
            q(i,j)=0.0;

	for (size_t i=0; i<m; i++)
		q(i,i)=1.0;
	for (size_t i=0; i<m; i++)
		for (size_t j=0; j<n; j++)
			r(i,j)=0.0;
	float_t ti;
	for (int i=d;i>0 && i--;)
	{
		vector_f h(m-i);
		for (size_t j=i;j<m;j++)
			h[j-i]=qr(j,i);
		matrix_f mm(m-i,m-i);
		for (size_t ii=i;ii<m;ii++)
			for (size_t jj=i; jj<m;jj++)
				mm(ii-i,jj-i)=q(ii,jj);
		ti=tau[i];
		householder_hm(ti,h,mm);
		for (size_t ii=i;ii<m;ii++)
			for (size_t jj=i; jj<m;jj++)
				q(ii,jj)=mm(ii-i,jj-i);
	}
	for (size_t i=0;i<m;i++)
		for (size_t j=i; j<n; j++)
			r(i,j)=qr(i,j);
}

void TFJSmoothing::qr(const matrix_f& a,
					matrix_f& q, matrix_f& r,
					bool ecnomic)
{
	matrix_f b=a;
	vector_f tau;
	qr_decomp(b,tau);
	qr_unpack(b,tau,q,r);

	if (ecnomic)
	{
		size_t m=a.size1(),
				n=a.size2();
		if (m<=n) return;
		matrix_f qq=q,rr=r;
		q.resize(m,n);
		r.resize(n,n);
		for (size_t i=0; i<m; i++)
			for (size_t j=0; j<n; j++)
				q(i,j)=qq(i,j);
		for (size_t i=0; i<n; i++)
			for (size_t j=0; j<n; j++)
				r(i,j)=rr(i,j);
	}
}

void TFJSmoothing::add_ends(const vector_f &c,
							vector_f &csm)
{
    using namespace std;
	size_t f=n_left+n_right+1;
	size_t hf=n_left;
	size_t od=n_order+1;
	matrix_f v(f,od);
	for (size_t i=0; i<f; i++)
		for (size_t j=0; j<od; j++)
			v(i,j)=1.0;
	vector_f t(f);
	for (int i=0; i<f; i++)
		t[i]=(float_t)i-hf;
	for (size_t i=0; i<od-1; i++)
		for (size_t j=0; j<f; j++)
			v(j,i+1)=pow(t[j],(float_t)(i+1));
	qr(v,_q,_r,true);

	vector_f bb0=begin_pts(c);
	for (size_t i=0; i<bb0.size(); i++)
		csm[i]=bb0[i];
	vector_f bb1=end_pts(c);
	for (size_t i=0; i<bb1.size(); i++)
		csm[csm.size()-bb1.size()+i]=bb1[i];
}

vector_f TFJSmoothing::begin_pts(const vector_f &c)
{
	size_t f=n_left+n_right+1;
	size_t hf=n_left;
	size_t od=n_order+1;

	matrix_f q=_q;
	// f=m; od=n;
	size_t m=q.size1(),
		   n=q.size2();

	// q(1:hf)*q'
	matrix_f p(n,m);
	for (size_t i=0; i<n; i++)
		for (size_t j=0; j<m; j++)
			p(i,j)=q(j,i);
	matrix_f qq(hf,f);
	for (size_t i=0; i<hf; i++)
		for (size_t j=0; j<f; j++)
		{
			float_t sm=0.0;
			for (size_t k=0; k<od; k++)
				sm+=q(i,k)*p(k,j);
			qq(i,j)=sm;
		}

	// *y(1:hf)
	vector_f d(hf);
	for (size_t i=0; i<hf; i++)
	{
		float_t sm=0.0;
		for (size_t j=0; j<f; j++)
			sm+=qq(i,j)*c[j];
		d[i]=sm;
	}
	return d;
}

vector_f TFJSmoothing::end_pts(const vector_f &c)
{
	size_t f=n_left+n_right+1;
	size_t hf=n_left;
	size_t od=n_order+1;

	matrix_f q=_q;
	// f=m; od=n;
	size_t m=q.size1(),
		   n=q.size2();

	// q(hf+2:end,:)*q'
	matrix_f p(n,m);
	for (size_t i=0; i<n; i++)
		for (size_t j=0; j<m; j++)
			p(i,j)=q(j,i);
	matrix_f qq(hf,f);
	for (size_t i=0; i<hf; i++)
		for (size_t j=0; j<f; j++)
		{
			float_t sm=0.0;
			for (size_t k=0; k<od; k++)
				sm+=q(hf+1+i,k)*p(k,j);
			qq(i,j)=sm;
		}

	// *y(1:hf)
	vector_f d(hf);
	for (size_t i=0; i<hf; i++)
	{
		float_t sm=0.0;
		for (size_t j=0; j<f; j++)
			sm+=qq(i,j)*c[c.size()-f+j];
		d[i]=sm;
	}
	return d;
}

}

