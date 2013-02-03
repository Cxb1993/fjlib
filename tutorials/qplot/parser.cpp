#include "fjlib_cio.h"
#include "fjlib_io.h"
#include "fjlib_vecmat.h"
#include "fjlib_string.h"
#include <iostream>
#include <fstream>

using namespace fjlib;

class TQPMatrix_Parser {
public:
	typedef std::vector<str_t> str_list_t;
	str_list_t	comments;
	matrix_f	mat;
};

inline
bool is_num(char a)
{
	return ((a>='0') && (a<='9'));
}

inline 
bool is_command(char a)
{
	return ((a=='!') || (a=='#') || (a=='['));
}


#include "fjlib_vecmat_print.h"
std::istream& operator>>(std::istream& in, TQPMatrix_Parser& p)
{
	// reset matrix
	p.comments.resize(0);

	char c[256];
	std::vector<std::vector<float> > data;
	size_t rows=0,cols=0;
	float v;
	char ch;
	// read all data lines
	while (!in.eof())
	{
		in.getline(c,256);
		if (in.eof()) break;
//		cout << c << endl;
		// detect comment or command line
		if (is_command(c[0])) {
			p.comments.push_back(c);
			continue;
		}

		data.resize(data.size()+1);
		std::vector<float> &row=data.back();
	
		std::istringstream ins(c);
		while (!ins.eof())
		{
			ch=ins.peek();
			// make sure it's number
			if (!is_num(ch)) {	
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
	
	return in;
}

int main()
{
	TQPMatrix_Parser a;
	quickload("in.txt",a);
//	TQPMatrix_Parser::str_list_t& cm=a.comments;
//	for_each(cm.begin(),cm.end(),print<str_t>);
	print_list("comments",a.comments);
	print_mat("matrix",a.mat);
}
