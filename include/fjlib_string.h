#ifndef FJLibStringH
#define FJLibStringH

#include <cstring> 
#include <sstream>
#include <iostream>

namespace fjlib {

typedef std::string str_t;

template< class type>
inline std::string to_string(const type & value)
{
    std::ostringstream streamOut;
    streamOut << value;
    return streamOut.str();
}

}

#endif

