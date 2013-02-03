#ifndef TFJLIBGroupH
#define TFJLIBGroupH

#include <vector>

namespace fjlib {

/*!
// Wrapper class to represent xyz group data. 
// It stores the data
*/
template <class T>
class TFJGroup_XYZ {
public:
	typedef T				value_type;
	typedef std::vector<T>	vec_type;
protected:
	/// Vector of type T, storage
	vec_type	vec;
public:
	TFJGroup_XYZ() {}
	/// Default constructor, with _dim dimension
	TFJGroup_XYZ(size_t _dim) { vec.resize(_dim); } 
	/// Returns the dimension of the group
	size_t		dim() const { return vec.size(); }
	/// Reset the dimension of the group
	void		redim(size_t _dim) { vec.resize(_dim); }

	/// Access the raw data
	vec_type&	vector() { return vec; }
	/// Return the data at ith node
	inline
	T&			at(size_t i) { return vec[i]; }
	/// Return the data at dimension x
	inline 
	T&			x() { return vec[0]; }
	/// Return the data at dimension y
	inline
	T&			y() { return vec[1]; }
	/// Return the data at dimension z
	inline
	T&			z() { return vec[2]; }

	/// Return the data at ith node
	inline const
	T&			at(size_t i) const { return vec[i]; }
	/// Return the data at dimension x
	inline const
	T&			x() const { return vec[0]; }
	/// Return the data at dimension y
	inline const
	T&			y()  const { return vec[1]; }
	/// Return the data at dimension z
	inline const
	T&			z() const { return vec[2]; }

};

template <class T>
std::ostream& operator<<(std::ostream& out, TFJGroup_XYZ<T>& c)
{
	using std::endl;
	out << c.dim() << endl; 
	for (size_t i=0; i<c.dim(); i++)
		out << c.at(i) << endl;
	return out;
}


} // end of namespace

#endif
