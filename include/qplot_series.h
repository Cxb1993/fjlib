#ifndef QPLOT_SERIES_H
#define QPLOT_SERIES_H

#include "qplot_matrix.h"

namespace qplot {

enum TQPSeriesType {
	qpsLine=0,qpsColormap
};

// created 3/8/06
struct TQPSeries {
	typedef TQPSeriesType series_type;
	TQPMatrix	data;
	bool		col_major;
	series_type	type;
	std::string	filename;
	
	TQPSeries(): col_major(true), type(qpsLine) {}
	// extracting command from comments and set rules
	void		ProcessComments() {
		if (data.comments.size()==0) return;	// no comments
		// need to be implemented later
	}
	// depending on the data structure, set default rules
	void		SetDefaultRules() {
		using namespace fjlib;
		matrix_f& m=data.mat;
		if ((m.size1()>2) && (m.size2()>2)) {
			type=qpsColormap;
			return;
		}
		type=qpsLine;
		col_major=(m.size1()>m.size2());
	}
};


}	// end of namespace

#endif
