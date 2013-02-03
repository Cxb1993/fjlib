#include "fjapp_MInterface.h"

namespace fjlib {

void TFJBubble_SInterface::import_pos(const vector_f& x,
								const vector_f& y)
{
	int n=x.size();
	_data.resize(MI_DATA_ROWS,n);
	for (size_t i=0; i<n; i++)
	{
		_data(drX,i)=x[i];
		_data(drY,i)=y[i];
	}
}

void TFJBubble_SInterface::export_pos(vector_f& x, vector_f& y)
{
	int n=_data.size2();
	x.resize(n); y.resize(n);
	for (size_t i=0; i<n; i++)
	{
		x[i]=_data(drX,i);
		y[i]=_data(drY,i);
	}
}

void TFJBubble_SInterface::import_row(MIDataRowType r, const vector_f& v)
{
	for (size_t i=0; i<size(); i++)
		_data(r,i)=v[i];
}

void TFJBubble_SInterface::export_row(MIDataRowType r, vector_f& v)
{
	v.resize(size());
	for (size_t i=0; i<size(); i++)
		v[i]=_data(r,i);
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TFJBubble_MInterface::export_pos(vector_f& x, vector_f& y)
{
	x.resize(0); y.resize(0);
	for (iter_type it=_mlist.begin(); it!=_mlist.end(); it++)
	{
		matrix_f &m=it->data();
		for (size_t i=0; i<m.size2(); i++)
		{
			push_back(x,m(drX,i));
			push_back(y,m(drY,i));
		}
	}
}

void TFJBubble_MInterface::export_row(MIDataRowType r, vector_f& v)
{
	v.resize(0);
	for (iter_type it=_mlist.begin(); it!=_mlist.end(); it++)
	{
		matrix_f &m=it->data();
		for (size_t i=0; i<m.size2(); i++)
			push_back(v,m(r,i));
	}
}

}	// end of namespace
