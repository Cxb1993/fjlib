#ifndef TFJLibPrintH
#define TFJLibPrintH

#include <iomanip>
#include "fjlib_string.h"
#include <fstream>

namespace fjlib {

template <class T>
class TFJPrint {
protected:
	T*				v_ptr;
	int				options;
	bool			head;
public:
	TFJPrint(T& v, bool w_head, int op=0): 
					v_ptr(&v), head(w_head), options(op) {}
	void			set(int op) { options=op; }
	virtual std::ostream&	
					ToString(std::ostream& os)
	{
		if (v_ptr==NULL) 
			os << "Null Pointer";
		else
			os << *v_ptr << std::endl;
		return os;
	}
	virtual std::ofstream&
					ToFile(std::ofstream& of)
	{
		if (v_ptr!=NULL) 
			of << *v_ptr << std::endl;
		return of;
	}
	virtual std::ifstream&
					LoadFile(std::ifstream& ifs)
	{
		if (v_ptr!=NULL)
			ifs >> *v_ptr; 
		return ifs;
	}
};

}	// end of namespaec

#endif
