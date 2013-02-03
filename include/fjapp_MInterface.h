//---------------------------------------------------------------------------
// Global Interface for fjISibubble
// It faciliates the handling of multiple interfaces


#ifndef FJAPP_MULTIPLE_INTERFACE
#define FJAPP_MULTIPLE_INTERFACE
//---------------------------------------------------------------------------
#include "fjlib.h"
#include <list>

#ifndef _MULTIPLE_INTERFACE
#define _MULTIPLE_INTERFACE
#endif

namespace fjlib {

//#define MI_DATA_ROWS 2		// clean
#define MI_DATA_ROWS 3			// surfactants

enum MIDataRowType {
	drX=0,drY,drGama
//	drXp1,drYp1,
//	drXp2,drYp2,
//	drUn,drUt,
//	drCurvature,
//	drGama			// index=9
};

/*!
// 10.22 created for multiple interface
// 10.23 arbitary  data import export supportedi
// 12.25 bug fixed due to resize of boost ublas matrix in linux
!*/
class TFJBubble_SInterface {
protected:
	void		init() { _data.resize(MI_DATA_ROWS,0); }
	matrix_f	_data;
public:
	TFJBubble_SInterface() { init(); }
	matrix_f&	data() { return _data; }
	const matrix_f&
				data() const { return _data; }

	// import, will resize automatically
	void		import_pos(const vector_f& x,
								const vector_f& y);
	void		export_pos(vector_f& x, vector_f& y);

	// import row, won't resize
	void		import_row(MIDataRowType r, const vector_f& v);
	void		export_row(MIDataRowType r, vector_f& v);

/*
	void		split(TFJBubble_SInterface& s1,
						TFJBubble_SInterface& s2, 
						size_t s1_count, bool joint_copy=false)
	{  
		fjlib::split(_data,s1._data,s2._data,
					s1_count,rcCol,joint_copy); 
	}
*/
	void		split(TFJBubble_SInterface& s,
						size_t s_count, bool joint_copy=false)
	{
		fjlib::split(_data,s._data,s_count,rcCol,joint_copy);
	}
	void		merge(const TFJBubble_SInterface& s, int merge_rule)
	{
		switch (merge_rule) {
		case 0:		// just regular merge
			fjlib::merge(_data,s._data,rcCol);
			break;	
		case 1:		// delete both the jointing node
			int n=_data.size2();
			fjlib::merge(_data,s._data,rcCol);
			fjlib::erase(_data,n-1,n,rcCol);	// not efficient!!
			break;
		}
	}
	void		erase(size_t pos1, size_t pos2)
	{
		fjlib::erase(_data,pos1,pos2,rcCol);
	}
	inline
	size_t		size() { return _data.size2(); }
	inline
	void		resize(size_t n) 
	{ 
		return fjlib::preserve_resize(_data,MI_DATA_ROWS,n); 
	}
};

class TFJBubble_MInterface {	
public:
	typedef TFJBubble_SInterface si_type;
	typedef std::list<si_type> list_type;
	typedef list_type::iterator iter_type;
protected:
	list_type	_mlist;
public:
	TFJBubble_MInterface() {}
	// get interface number
	inline
	size_t			size() { return _mlist.size(); }
	void			clear() { _mlist.clear(); }

	inline
	iter_type		append() 
	{ 
		si_type a; return _mlist.insert(_mlist.end(),a);	
	}
	inline
	iter_type		insert(iter_type it)
	{ si_type a; return _mlist.insert(it,a); }

	inline
	list_type&		mlist() { return _mlist; }

	// assemble the coordinates for output
	void			export_pos(vector_f& x, vector_f& y);
	// assemble other properties
	void			export_row(MIDataRowType r, vector_f& v);
};

}

#endif
