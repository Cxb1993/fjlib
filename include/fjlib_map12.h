//---------------------------------------------------------------------------

#ifndef FJMapH
#define FJMapH

#include <functional>
#include <map>

namespace fjlib {

template <class _Key, class _Tp1, class _Tp2,
class _Compare=std::less<_Key>,
class _Allocator=std::allocator<std::pair<const _Key,
	std::pair<_Tp1,_Tp2> > > >
struct mapp {
        typedef _Key                            key_type;
        typedef _Tp1                            data1_type;
        typedef _Tp2                            data2_type;
		typedef std::pair<_Tp1,_Tp2>            data_type;
		typedef std::pair<_Tp1,_Tp2>            mapped_type;
		typedef std::pair<const _Key,data_type> value_type;
        typedef _Compare                        key_compare;
		typedef std::map<_Key,data_type,_Compare,_Allocator>	
												map_type;
        typedef typename map_type::iterator     iterator;
        typedef typename map_type::const_iterator
												const_iterator;
};

}  // end of namespace



//---------------------------------------------------------------------------
#endif
