#ifndef FJFEMStorageH
#define FJFEMStorageH

// 11.16 seperate storage from element and build
// 1.3	 use list and map to access node and element list

#include "fjlib_fem_element.h"
#include <list>
#include <map>

namespace fjlib {

// TFJFEMStorage allows the transformation from
// the raw data to the linked list of the nodes and elements,
// including nodes location and boundary coefficient.
// It's supposed to be fed into TFJFEMBuilder
class TFJFEMStorage_Base {
public:
	virtual	
	void			build()=0;			// build elements list
	virtual
	size_t			nodes_count()=0;	// nodes count
	virtual
	size_t			elements_count()=0;	// elements count
};

/*
// Boundary usually requires the data of geometry
// therefore it's defined here instead of "fjlib_fem_element.h"
struct TFJFEMBoundaryData: public TFJFEMElementData {
	int				element_id;
	int				bc_type;
	TFJFEMStorage_Base*
					storage;
};
*/

// Node Structure
// includes all the data about one node
// info of the element can be only accessed after assembled
struct TFJFEMNodeStore {
	int			node_id;			// index of node 
									// don't have to be in order
	int			matrix_id;			// reserved for matrix assembly
	void*		element;			// where it belong to
	float_t		pos[3];				// coordinates
	int			bc_type;			// boundary type
									// -1, no bc
	vector_f	bc_coefs;			// boundary coefficients
									// reserved
};

// Element Structure
// includes all the data about one elemeent
// linked with all the nodes it's composed of
struct TFJFEMElementStore {
	int			element_id;			// index of element
	TFJFEMElementSpec*
				spec;				// element specification
	std::vector<TFJFEMNodeStore*>
				nodes;				// nodes link list
	std::vector<float_t>
				bc_types;			// boundary types
									// zero size, no bc
	vector_f	bc_coefs;			// reserved
};

class TFJFEMStorage: public TFJFEMStorage_Base {
public:
	typedef TFJFEMNodeStore node_type;
	typedef TFJFEMElementStore element_type;
	typedef std::list<node_type> node_list_type;
	typedef std::list<element_type> element_list_type;
	typedef TFJFEMElementSpec spec_type;
	typedef std::vector<spec_type> spec_vec_type;
public:
	node_list_type	nodes;
	element_list_type
					elements;
protected:
	spec_vec_type*	specs;
	matrix_f*		nlist;					// node pos list
	matrix_n*		elist;					// element nodes list
	matrix_f*		node_bc;				// boundary type 1 list
	matrix_f*		el_bc;					// boundary type 2 list
public:
	typedef std::map<int,void*> map_type;
	typedef std::pair<int,void*> pair_type;
protected:
	map_type		node_map,				// mapping from
					el_map;					// node_id to 
											// node pointer
public:
	node_type*		get_node(int node_id)
	{ return (node_type*)node_map[node_id]; }
	element_type*	get_element(int el_id)
	{ return (element_type*)el_map[el_id]; }
public:
	void			set(spec_vec_type* sp, matrix_f* nl, 
						matrix_n* el, matrix_f* nbc,
						matrix_f* elbc)
	{ specs=sp; nlist=nl; elist=el; node_bc=nbc; el_bc=elbc; }
	void			build();			// build elements list
	size_t			nodes_count() { return nodes.size(); }
	size_t			elements_count() { return elements.size(); }
/*
	void			get_nodes(int el_id,
							vector_n* nodes_storage_id);
*/
public:
	// utilities, helper function
	void			get_sizes(element_type* el, vector_f* sz);
	void			get_matrix_ids(element_type* el, 
								vector_n* ids);
	void			get_nodes_pos(element_type* el, matrix_f* pos);
protected:
	matrix_f*		solution;		// solution or property
public:
	void			set_solution(matrix_f* sl)
	{ solution=sl; }
	matrix_f*		get_solution() 
	{ return solution; }
	bool			has_solution() { return (solution!=NULL); }
	std::ofstream&	ToFile(std::ofstream& of);
public:
	void			get_nodes_value(vector_n* ids, matrix_f* v);
public:
	TFJFEMStorage():solution(NULL) {}
};

}	// end of namespace

#endif