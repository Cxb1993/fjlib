#ifndef QPLOT_MATRIX_PARSER_H
#define QPLOT_MATRIX_PARSER_H

#include "fjlib_string.h"
#define USE_STL_MATRIX
#include "fjlib_vecmat.h"

namespace qplot {

// created 3/6/06
// 3/22/06 support two format, auto or manual
// auto, just row x col data
// manual, one line indicate total, ex
//		for	vector, 
//			10
//			1 2 3 4 5 6 7 8 9 10
//			5
//			1 2 3 4 5
//		for matrix,
//			2 3			or 		2 3
//			1 2 3 				1 2 
//			4 5 6				3 4 
//								5 6
struct TQPMatrix {
	typedef std::vector<fjlib::str_t> str_list_t;
	str_list_t		comments;
	fjlib::matrix_f	mat;
};

inline
bool _is_num(char a)
{
	return (((a>='0') && (a<='9')) || (a=='-') || (a=='.'));
}

inline 
bool _is_command(char a)
{
	return ((a=='!') || (a=='#') || (a=='['));
}

inline
bool _is_integer(fjlib::float_t v)
{
	int n=(int)v;
	return (((v-n)==0) && (v>0));
}

void _throw_error()
{
	throw "format not supported!";
}

void _process_mat(std::vector<std::vector<
							fjlib::float_t> >& raw,
					fjlib::matrix_f& m)
{
	bool done=false;
	size_t i,j,k,rows,cols,ib;
	
	// figure out the structure
	fjlib::matrix_sz index(raw.size(),3); 	// 0 rows
										  	// 1 cols
											// 2 next physical row
	for (i=0; i<raw.size(); i++) index(i,0)=0;
	i=0;
	while ((!done) && (i<raw.size())) {
		std::vector<fjlib::float_t> &row=raw[i];
		// check expect type
		bool type=0;					// auto
		if (row.size()==1) 
			if (_is_integer(row[0])) {
				type=1;					// vector
				rows=1; cols=(int)row[0];
			}
		if (row.size()==2)
			if (_is_integer(row[0]) && _is_integer(row[1])) {
				type=2;					// matrix
				rows=(int)row[0]; cols=(int)row[1];
				done=true;
			}
		if (type==0) {
			rows=1; cols=row.size();
		}
		size_t tot;
		ib=i;
		switch (type) {
			case 0:
				i++;
				break;
			default: // 1 or 2
				index(i,0)=0;
				i++; if (i==raw.size()) _throw_error();
				ib=i; // start point
				j=0; tot=rows*cols;
				// search tot number of numbers
				while ((j<tot) && (i<raw.size())) {
					j+=raw[i].size();
					index(i,0)=0;
					i++;
				}
				if (j<tot) _throw_error();
				break;
		}
		index(ib,0)=rows;
		index(ib,1)=cols;
		index(ib,2)=i;
	}
	// get max rows and cols
	i=0; rows=0; cols=0;
	while (i<raw.size()) {
		rows+=index(i,0);
		if (index(i,0)>0) {
			if (index(i,1)>cols) cols=index(i,1);
		}
		i++;
	}
	// make final matrix 
	m.resize(rows,cols);
	i=0; 	// current raw row 
	k=0;	// mat row
	while (i<raw.size()) {
		if (index(i,0)<1) { i++; continue; }
		
		ib=i; j=0; done=false;
		for (size_t ii=0; ii<index(ib,0); ii++) {
			if (done) break;
			for (size_t jj=0; jj<index(ib,1); jj++) {
				m(k,jj)=raw[i][j];
				j++; 
				if (j==raw[i].size()) { 
					j=0; i++;
					if (i==index(ib,2)) { done=true; break; }
				}
			}
			k++;
		}
		i=index(ib,2);	// go to next avariable row
	}
}
					

std::istream& operator>>(std::istream& in, TQPMatrix& p)
{
	// reset matrix
	p.comments.resize(0);

	std::vector<std::vector<fjlib::float_t> > data;
	size_t rows=0,cols=0;
	fjlib::float_t v;
	fjlib::str_t c;
	char ch;
	// read all data lines
	while (!in.eof())
	{
		c="";
		// get a line
		do {
			ch=in.get();
			if (ch<0) break;	// end of file
			if (ch!='\n') 
				c.append(&ch,1);
			else break;
		} while (!in.eof());
		if (c=="") continue;
//		in.getline(c,256);
//		if (c=='\0') continue;
//		cout << c << endl;
		// detect comment or command line
		if (_is_command(c[0])) {
			p.comments.push_back(c);
			continue;
		}

		data.resize(data.size()+1);
		std::vector<fjlib::float_t> &row=data.back();
	
		std::istringstream ins(c);
		while (!ins.eof())
		{
			ch=ins.peek();
			// make sure it's number
			if (!_is_num(ch)) {	
				ins.ignore(1);
				continue;
			}
			ins >> v;
			row.push_back(v);
		}
		if (row.size()==0) continue;
		if (row.size()>cols) cols=row.size();
		rows++;
	}
	// process it
	_process_mat(data,p.mat);
/*	
	p.mat.resize(rows,cols);
	// copy data to a matrix
	size_t i=0;
	while (i<rows) {
		size_t j=0;
		while (j<data[i].size()) {
			p.mat(i,j)=data[i][j];
			j++;
		}
		i++;
	}
*/	
	return in;
}

}	// end of namespace

#endif
